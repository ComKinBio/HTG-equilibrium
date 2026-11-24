"""Pin the Julia version and every registered Julia package in src/htg/juliapkg.json.

juliapkg (used by pyclapeyron/juliacall) otherwise picks the newest Julia and
the newest package versions allowed by the bounds in the installed packages,
so a fresh environment can resolve different versions than the one the
reference results were checked with. Standard libraries (LinearAlgebra,
OpenBLAS_jll, ...) come with the Julia version, so pinning Julia pins them.

Usage, after checking that the tests pass with the current environment:

    uv run python scripts/pin_julia.py [path/to/Manifest.toml]

The default manifest is the one juliapkg resolved for this virtual environment
(.venv/julia_env/Manifest.toml). The script rewrites src/htg/juliapkg.json
with exact ("=") versions taken from that manifest, keeping the package's own
direct dependencies first.

To update the Julia environment: loosen or remove the pins in
src/htg/juliapkg.json, run any calculation (juliapkg resolves again), run the
tests, then run this script to pin the new versions.
"""
import json
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEPS_FILE = ROOT / 'src' / 'htg' / 'juliapkg.json'
DEFAULT_MANIFEST = ROOT / '.venv' / 'julia_env' / 'Manifest.toml'
# Pinned exactly by juliacall, whose version uv.lock fixes.
OWNED_ELSEWHERE = {'PythonCall'}


def pins(manifest_path):
    manifest = tomllib.loads(Path(manifest_path).read_text())
    packages = {}
    for name, entries in manifest['deps'].items():
        (entry,) = entries
        # Registered packages have a git tree hash; standard libraries do not
        # and are fixed by the Julia version.
        if 'git-tree-sha1' not in entry or name in OWNED_ELSEWHERE:
            continue
        # Compat specifiers have no build numbers (JLL "1.2.3+0" -> "=1.2.3").
        version = entry['version'].split('+')[0]
        packages[name] = {'uuid': entry['uuid'], 'version': '=' + version}
    return manifest['julia_version'], packages


def main():
    manifest_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_MANIFEST
    julia, packages = pins(manifest_path)
    current = json.loads(DEPS_FILE.read_text())
    direct = list(current.get('direct', current['packages']))
    ordered = {n: packages[n] for n in direct if n in packages}
    ordered |= {n: packages[n] for n in sorted(packages, key=str.lower) if n not in ordered}
    data = {
        'julia': '=' + julia,
        'direct': direct,
        'packages': ordered,
    }
    DEPS_FILE.write_text(json.dumps(data, indent=4) + '\n')
    print(f'Pinned Julia {julia} and {len(ordered)} packages in {DEPS_FILE.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
