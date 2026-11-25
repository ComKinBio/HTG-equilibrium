"""Write test/data/state_python.json: what the Python app's State and status line hold.

    uv run python julia/HTGApp/test/dump_state.py

- the species lists of State.use_model for the shipped data (no user species):
  phase_options (gas, aqueous, organic_liquid, solid) and feed_species
- worker_summary for a set of worker states and queue lengths
"""
import json
from pathlib import Path

from htg.gui.common import worker_summary
from htg.model import Model

HERE = Path(__file__).resolve().parent
OUT = HERE / 'data' / 'state_python.json'

STATUSES = [
    ([], 0),
    (['starting'] * 4, 0),
    (['warming', 'idle', 'idle', 'idle'], 0),
    (['busy', 'busy', 'idle'], 0),
    (['busy', 'busy', 'busy'], 7),
    (['crashed', 'idle', 'stopping', 'busy', 'warming', 'starting'], 1),
    (['idle'], 2),
    ([], 3),
]


def main():
    model = Model.load(user_data=None)
    catalogue = model.species()
    molecular = catalogue[catalogue.Phase == 'molecular'].Species.tolist()
    ions = catalogue[catalogue.Phase == 'aqueous ion'].Species.tolist()
    solids = catalogue[catalogue.Phase == 'solid'].Species.tolist()
    registry = model.registry
    summaries = []
    for states, queued in STATUSES:
        status = {'workers': [{'id': i + 1, 'state': s, 'message': '', 'case': '', 'point': None}
                              for i, s in enumerate(states)],
                  'queued_pieces': queued, 'target': len(states)}
        summaries.append({'states': states, 'queued_pieces': queued, 'text': worker_summary(status)})
    data = {
        'phase_options': {'gas': molecular, 'aqueous': molecular + ions, 'organic_liquid': molecular,
                          'solid': solids},
        'feed_species': sorted({*registry.all_molecular, *registry.solid}),
        'worker_summary': summaries,
    }
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(data, indent=0, ensure_ascii=False) + '\n', encoding='utf-8')
    print('Wrote', OUT)


if __name__ == '__main__':
    main()
