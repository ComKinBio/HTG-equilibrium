"""Command line interface.

    htg run CASE.toml [CASE2.toml ...] [--out DIR] [--workers N] [--no-plots]
    htg check CASE.toml
    htg species [--all] [--phase PHASE] [--search TEXT]
    htg new CASE.toml
    htg plot RESULTS_DIR
    htg gui [--workers N] [--port PORT] [--results DIR] [--cases DIR] [--no-browser]
"""
import argparse
import sys
import time
from pathlib import Path

import pandas as pd

TEMPLATE = Path(__file__).parent / 'templates' / 'case_template.toml'


def _model(args):
    from .model import Model
    return Model.load(args.data)


def _print_table(frame, columns=None):
    with pd.option_context('display.max_rows', None, 'display.max_columns', None,
                           'display.width', 200, 'display.max_colwidth', 80):
        print((frame[columns] if columns else frame).to_string(index=False))


def cmd_run(args):
    from .case import load_case
    cases = [load_case(path) for path in args.case]
    names = [c.name for c in cases]
    if len(set(names)) != len(names):
        raise ValueError(f'Case names must differ (they name the output folders): {names}')
    if args.out and len(cases) == 1:
        outputs = [Path(args.out)]
    else:
        outputs = [Path(args.out or 'results') / case.name for case in cases]
    model = _model(args)
    for case in cases:
        points = len(case.grid())
        print(f'Case {case.name}: {points} point(s), {len(case.chains())} independent piece(s), '
              f'solver {case.solver}')
        check = model.check(case)
        invalid = check[check.status != 'ready']
        if len(invalid):
            print('Preflight found problems at some points (they will probably fail):')
            _print_table(invalid, ['T_C', 'P_bar', 'error'])
    start = time.perf_counter()
    if args.workers > 1:
        from .runner import WorkerPool

        def on_point(run, T, P, outcome, start_label):
            text = 'FAILED: ' + outcome if isinstance(outcome, str) else 'solved'
            print(f'{run.name}: {T:g} °C, {P:g} bar: {text} (start: {start_label})', flush=True)
        print(f'Starting {args.workers} worker processes (each starts Julia; the first results '
              'take about a minute)')
        with WorkerPool(args.workers, model.load_args, warm_up=False, on_point=on_point) as pool:
            runs = [pool.submit(case, out) for case, out in zip(cases, outputs)]
            for run in runs:
                run.wait()
        sweeps = [run.sweep_result() for run in runs]
    else:
        sweeps = []
        for case, out in zip(cases, outputs):
            sweep = model.sweep(case)
            sweep.save(out)
            sweeps.append(sweep)
    print(f'Finished in {time.perf_counter() - start:.0f} s')
    ok = True
    for case, sweep, out in zip(cases, sweeps, outputs):
        summary = sweep.summary()
        if 'error' in summary:
            summary['error'] = summary.error.fillna('').str.slice(0, 90)
        print(f'\n{case.name}:')
        _print_table(summary, [c for c in ['T_C', 'P_bar', 'status', 'gibbs_J_per_mol_C',
                                           'gas_carbon_percent', 'aqueous_removed', 'error']
                               if c in summary])
        if not args.no_plots:
            _plots(out, case)
        print(f'Results written to {out}/')
        ok = ok and bool((summary.status == 'solved').all())
    return 0 if ok else 1


def _plots(directory, case=None):
    import matplotlib
    matplotlib.use('Agg')
    from .plotting import save_plots
    comparison, quantity, label = None, None, 'Reference'
    if case is not None and case.comparison:
        comparison = pd.read_csv(case.comparison['file'])
        quantity = case.comparison.get('quantity')
        label = case.comparison.get('label', label)
    for path in save_plots(directory, quantity, comparison, label):
        print(f'Plot: {path}')


def cmd_plot(args):
    from .case import load_case
    directory = Path(args.results)
    case_file = directory / 'case.toml'
    _plots(directory, load_case(case_file) if case_file.exists() else None)
    return 0


def cmd_check(args):
    from .case import load_case
    case = load_case(args.case)
    table = _model(args).check(case)
    _print_table(table)
    ok = (table.status == 'ready').all()
    print('All temperatures passed the preflight checks.' if ok else
          'Some temperatures are invalid; see the error column.')
    return 0 if ok else 1


def cmd_species(args):
    table = _model(args).species(include_unavailable=args.all, phase=args.phase, search=args.search)
    columns = ['Species', 'Formula', 'Phase'] + (['Ready', 'Reason'] if args.all else [])
    _print_table(table, columns)
    print(f'{len(table)} species')
    return 0


def cmd_new(args):
    path = Path(args.case)
    if path.exists():
        print(f'{path} already exists; choose another name.', file=sys.stderr)
        return 1
    text = TEMPLATE.read_text(encoding='utf-8').replace('name = "my_case"', f'name = "{path.stem}"')
    path.write_text(text, encoding='utf-8')
    print(f'Wrote {path}. Edit it, then: htg check {path}  and  htg run {path}')
    return 0


def cmd_gui(args):
    from .gui.app import run_app
    run_app(results_dir=args.results, cases_dir=args.cases, workers=args.workers, port=args.port,
            show=not args.no_browser, host=args.host, model_args={'data_dir': args.data} if args.data else None)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog='htg', description='Hydrothermal gasification equilibrium.')
    parser.add_argument('--data', help='database folder (default: the repository data/ folder)')
    sub = parser.add_subparsers(dest='command', required=True)

    p = sub.add_parser('run', help='solve case files and write tables and plots')
    p.add_argument('case', nargs='+')
    p.add_argument('--out', help='output folder (one case), or parent folder (several cases); '
                                 'default: results/<case name>')
    p.add_argument('--workers', type=int, default=1,
                   help='worker processes; independent points and continuation chains of all '
                        'cases are shared between them')
    p.add_argument('--no-plots', action='store_true')
    p.set_defaults(func=cmd_run)

    p = sub.add_parser('check', help='fast checks of a case file (no equation of state)')
    p.add_argument('case')
    p.set_defaults(func=cmd_check)

    p = sub.add_parser('species', help='list species that can be used in case files')
    p.add_argument('--all', action='store_true', help='include species with incomplete data')
    p.add_argument('--phase', help="e.g. molecular, organic, solid, ion")
    p.add_argument('--search', help='part of a name or formula')
    p.set_defaults(func=cmd_species)

    p = sub.add_parser('new', help='write a commented case-file template')
    p.add_argument('case')
    p.set_defaults(func=cmd_new)

    p = sub.add_parser('plot', help='re-plot saved results')
    p.add_argument('results')
    p.set_defaults(func=cmd_plot)

    p = sub.add_parser('gui', help='open the app in the web browser')
    p.add_argument('--workers', type=int, help='worker processes (default: from free cores and memory)')
    p.add_argument('--port', type=int, default=8080)
    p.add_argument('--host', default='127.0.0.1',
                   help='network address (default 127.0.0.1: only this computer; 0.0.0.0: everyone on the network)')
    p.add_argument('--results', default='results', help='folder for results (default: results)')
    p.add_argument('--cases', default='cases', help='folder with case files (default: cases)')
    p.add_argument('--no-browser', action='store_true', help='do not open a browser window')
    p.set_defaults(func=cmd_gui)

    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, FileNotFoundError) as exc:
        print(f'Error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
