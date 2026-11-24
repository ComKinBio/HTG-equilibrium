"""Write test/data/python_cli/: what the Python command line prints, for the Julia tests.

    uv run python julia/HTGEquilibrium/test/dump_cli.py

Runs `htg.cli.main(argv)` in this process (stdout, stderr and exit code
captured) from the repository root, for:

- usage errors, help, and errors reported as "Error: ..." (exit 2);
- `htg check` of every case file, the old case files, the template and a few
  variants of the template with problems (written to a temporary folder; their
  text is stored);
- `htg species` with --all, --phase, --search, and with the user species folder
  of test/data/user_species;
- `htg new` (the file it writes);
- pandas `to_string(index=False)` of a few hand-made tables (every dtype, NaN,
  None, floats that pandas writes in exponent notation, an empty table);
- `htg run` of a small case (three independent points of
  yakaboylu2013_fig7_8, one of them failing): its output and the folder it
  writes (run/, CSV files and case.toml).

Long outputs (`species --all`, ~400 kB) are stored as sha256, line count and
first lines only. Words in argv: {USER_SPECIES} is test/data/user_species,
{TMP} a temporary folder holding the variants and the run case.
"""
import contextlib
import hashlib
import io
import json
import os
import re
import shutil
import tempfile
from pathlib import Path

from htg.cli import main

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE / 'data' / 'python_cli'
TEMPLATE = REPO / 'src' / 'htg' / 'templates' / 'case_template.toml'
USER_SPECIES = HERE / 'data' / 'user_species'
FULL_TEXT_LIMIT = 20000
NO_USER = '--no-user-species'

USAGE = [[], ['-h'], ['--help'], ['foo'], ['run'], ['run', '-h'], ['check', '-h'], ['species', '-h'],
         ['new', '-h'], ['--data'], ['--x', 'species'], ['species', '--x'], ['species', 'extra'],
         ['species', '--phase'], ['species', '--search', '--phase', 'ion'], ['species', '--all=1'],
         ['run', 'x.toml', '--workers', 'abc'], ['run', 'x.toml', '--workers'], ['run', '--workers=2'],
         ['run', 'a.toml', '--no-plots', 'b.toml'], ['check'], ['check', 'a', 'b'], ['new'], ['new', 'a', 'b'],
         ['gui', '--port', 'x'], ['plot'], ['species', '--no-user', '--user-data'],
         ['species', '--', '--all'], ['check', '--', '-h']]
ERRORS = [['check', 'missing.toml', NO_USER], ['check', './cases/../missing.toml', NO_USER],
          ['run', 'missing.toml', NO_USER], ['--data', '/nonexistent', 'species', NO_USER],
          ['run', 'cases/manure_gas_min.toml', 'cases/manure_gas_min.toml', NO_USER],
          ['new', 'cases/manure_gas_min.toml'], ['new', 'no_such_folder/x.toml'],
          ['check', '{TMP}/bad_syntax.toml', NO_USER], ['check', '{TMP}/bad_value.toml', NO_USER]]
SPECIES = [['species', NO_USER], ['species', '--all', NO_USER], ['species', '--phase', 'ion', NO_USER],
           ['species', '--phase', 'solid', NO_USER], ['species', '--phase', 'molecular', '--all', NO_USER],
           ['species', '--phase', 'SUPCRT', '--all', NO_USER], ['species', '--search', 'meth', NO_USER],
           ['species', '--search', 'CO2', '--all', NO_USER], ['species', '--search', 'zzz', NO_USER],
           ['species', '--search', 'zzz', '--all', NO_USER], ['species', '--phase', 'organic', '--search', 'ene', NO_USER],
           ['species', '--ph', '^(?:ion|solid)$', '--sea', 'Na', NO_USER], ['species', '--phase=aq', '--search=H', NO_USER],
           ['species', '--user-data', '{USER_SPECIES}'], ['species', '--user-data', '{USER_SPECIES}', '--all'],
           ['species', '--user-data', '{USER_SPECIES}', '--search', 'my', '--all']]


NAN, INF = '__nan__', '__inf__'   # JSON has neither
TABLES = [  # columns -> values; None is Python's None
    {'T_C': [300.0, 250.5, 1000.0], 'P_bar': [240.0, 240.0, 1.0], 'n': [1, 22, -3], 'ok': [True, False, True]},
    {'small': [1.5, NAN, 1e-7], 'big': [-481234.567891234, -3328077.1, 12.0], 'mixed': [1e7, 2.0, 0.1 + 0.2]},
    {'large_short': [1234567.0, 1.0, NAN], 'large_long': [123456789.123, 1.0, 2.0], 'edge': [1e-6, 1.0, -2.5e-6]},
    {'zeros': [0.0, -0.0, 5.0], 'inf': [INF, 1.0, -1.25], 'neg': [-1.0, -0.5, -100.125]},
    {'bool_none': [True, None, False], 'str_none': ['a', None, 'long text\nwith newline\tand tab'],
     'all_none': [None, None, None]},
    {'object': ['x', 1, 2.5], 'object2': [True, 'y', -1.5], 'unicode': ['°C', 'Ä', 'x' * 30],
     'with space': ['', '', '']},
    {'gibbs_J_per_mol_C': [-3328077.123456, NAN, -3457880.5], 'aqueous_removed': [False, NAN, True],
     'error': ['', 'FAILED: x, "quoted"', '']},
]


def _tables():
    import pandas as pd
    out = []
    for spec in TABLES + [{}]:
        data = {c: [float('nan') if v == NAN else float('inf') if v == INF else v for v in values]
                for c, values in spec.items()}
        frame = pd.DataFrame(data)
        if not spec:
            frame = pd.DataFrame(columns=['Species', 'Formula', 'Phase'])
        with pd.option_context('display.max_rows', None, 'display.max_columns', None,
                               'display.width', 200, 'display.max_colwidth', 80):
            text = frame.to_string(index=False)
        out.append({'columns': spec, 'dtypes': [str(frame[c].dtype) for c in frame.columns], 'text': text})
    return out


def _variants():
    """Template variants with preflight problems: name -> text."""
    text = TEMPLATE.read_text(encoding='utf-8')
    temps = 'temperature_C = [300, 350, 400, 450]'
    assert temps in text and '"nitrogen"]' in text and 'cutoff_C = 390.0' in text and 'solid = ["carbon"]' in text
    return {
        'check_unknown': text.replace('"nitrogen"]', '"nitrogen", "unobtainium"]'),
        'check_no_cutoff': text.replace('cutoff_C = 390.0', 'cutoff_C = 600.0'),
        'check_hot': text.replace(temps, 'temperature_C = [300, 3500]'),
        'check_mineral': text.replace(temps, 'temperature_C = [300, 1200]').replace(
            'solid = ["carbon"]', 'solid = ["carbon", "Halite"]'),
        'check_points': re.sub(r'temperature_C = .*\npressure_bar = [^\n]*', 'points = [[300, 240], [350.5, 250], [400, 1]]',
                               text),
        'bad_syntax': text.replace('[phases]', '[phases'),
        'bad_value': text.replace('pressure_bar = 240.0', 'pressure_bar = -5.0'),
    }


def _run_case():
    """The `htg run` case: three independent points of yakaboylu2013_fig7_8 (389.6 °C fails)."""
    text = (REPO / 'cases' / 'yakaboylu2013_fig7_8.toml').read_text(encoding='utf-8')
    text, n = re.subn(r'temperature_C = \[[^\]]*\]', 'temperature_C = [300, 389.6, 450]', text)
    assert n == 1
    return text.replace('name = "yakaboylu2013_fig7_8"', 'name = "cli_run_small"')


def run(argv, cwd):
    out, err = io.StringIO(), io.StringIO()
    old = os.getcwd()
    os.chdir(cwd)
    try:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            try:
                code = main(argv)
            except SystemExit as exc:
                code = exc.code
    finally:
        os.chdir(old)
    entry = {'argv': argv, 'code': code, 'stderr': err.getvalue()}
    text = out.getvalue()
    if len(text) > FULL_TEXT_LIMIT:
        entry['stdout_sha256'] = hashlib.sha256(text.encode('utf-8')).hexdigest()
        entry['stdout_lines'] = text.count('\n')
        entry['stdout_head'] = ''.join(text.splitlines(keepends=True)[:5])
    else:
        entry['stdout'] = text
    return entry


def main_dump():
    tmp = Path(tempfile.mkdtemp(prefix='htg_dump_cli_'))
    variants = _variants()
    for name, text in variants.items():
        (tmp / f'{name}.toml').write_text(text, encoding='utf-8')

    def expand(argv):
        return [a.replace('{TMP}', str(tmp)).replace('{USER_SPECIES}', str(USER_SPECIES)) for a in argv]

    def record(argv, cwd=REPO):
        entry = run(expand(argv), cwd)
        entry['argv'] = argv
        for key in ('stdout', 'stderr'):
            if key in entry:
                entry[key] = entry[key].replace(str(tmp), '{TMP}')
        return entry

    checks = ([['check', f'cases/{p.name}', NO_USER] for p in sorted((REPO / 'cases').glob('*.toml'))]
              + [['check', f'tests/data/legacy_cases/{p.name}', NO_USER]
                 for p in sorted((REPO / 'tests' / 'data' / 'legacy_cases').glob('*.toml'))]
              + [['check', 'src/htg/templates/case_template.toml', NO_USER]]
              + [['check', f'{{TMP}}/{name}.toml', NO_USER] for name in variants if name.startswith('check_')])
    commands = [record(argv) for argv in USAGE + ERRORS + checks + SPECIES]

    new_dir = tmp / 'new'
    new_dir.mkdir()
    new = run(['new', 'my_new_case.toml'], new_dir)
    new['file'] = (new_dir / 'my_new_case.toml').read_text(encoding='utf-8')
    (new_dir / 'sub.dir').mkdir()
    new_other = run(['new', 'sub.dir//other'], new_dir)
    new_other['file'] = (new_dir / 'sub.dir' / 'other').read_text(encoding='utf-8')

    run_dir = tmp / 'run'
    (run_dir / 'comparison').mkdir(parents=True)
    shutil.copy(REPO / 'cases' / 'comparison' / 'yakaboylu2013_fig7_8.csv', run_dir / 'comparison')
    case_text = _run_case()
    (run_dir / 'cli_run_small.toml').write_text(case_text, encoding='utf-8')
    run_argv = ['run', 'cli_run_small.toml', '--out', 'out', '--no-plots', NO_USER]
    run_entry = run(run_argv, run_dir)
    run_entry['stderr'] = ''   # juliacall's set-up messages

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    shutil.copytree(run_dir / 'out', OUT / 'run')
    data = {'tables': _tables(), 'commands': commands, 'variants': variants, 'new': [new, new_other],
            'run': {'case': case_text, 'argv': run_argv, **run_entry}}
    (OUT / 'cli_python.json').write_text(json.dumps(data, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
    shutil.rmtree(tmp)
    size = sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
    print(f'{OUT}: {len(commands)} commands, run of {run_argv[1]}, {size / 1000:.0f} kB')


if __name__ == '__main__':
    main_dump()
