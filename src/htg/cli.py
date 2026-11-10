"""Command line interface.

    htg run CASE.toml [--out DIR] [--workers N] [--no-plots]
    htg check CASE.toml
    htg species [--all] [--phase PHASE] [--search TEXT]
    htg new CASE.toml
    htg plot RESULTS_DIR
"""
import argparse
import shutil
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
    case = load_case(args.case)
    out = Path(args.out) if args.out else Path('results') / case.name
    print(f'Case {case.name}: {len(case.temperatures())} temperature(s), solver {case.solver}')
    model = _model(args)
    check = model.check(case)
    invalid = check[check.status != 'ready']
    if len(invalid):
        print('Preflight found problems at some temperatures (they will probably fail):')
        _print_table(invalid, ['T_C', 'error'])
    start = time.perf_counter()
    sweep = model.sweep(case, workers=args.workers)
    print(f'Finished in {time.perf_counter() - start:.0f} s')
    sweep.save(out)
    saved = case.replace()
    if case.comparison:
        shutil.copyfile(case.comparison['file'], out / 'comparison.csv')
        saved.comparison = {**case.comparison, 'file': 'comparison.csv'}
    saved.save(out / 'case.toml')
    summary = sweep.summary()
    if 'error' in summary:
        summary['error'] = summary.error.fillna('').str.slice(0, 90)
    _print_table(summary, [c for c in ['T_C', 'status', 'gibbs_J_per_mol_C', 'gas_carbon_percent',
                                       'aqueous_removed', 'error'] if c in summary])
    if not args.no_plots:
        _plots(out, case)
    print(f'Results written to {out}/')
    return 0 if (summary.status == 'solved').all() else 1


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


def main(argv=None):
    parser = argparse.ArgumentParser(prog='htg', description='Hydrothermal gasification equilibrium.')
    parser.add_argument('--data', help='database folder (default: the repository data/ folder)')
    sub = parser.add_subparsers(dest='command', required=True)

    p = sub.add_parser('run', help='solve a case file and write tables and plots')
    p.add_argument('case')
    p.add_argument('--out', help='output folder (default: results/<case name>)')
    p.add_argument('--workers', type=int, default=1,
                   help='parallel processes for sweeps with continuation = false')
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

    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, FileNotFoundError) as exc:
        print(f'Error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
