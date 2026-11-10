"""The worker pool's scheduling, with stand-in workers (no Julia), plus one real run."""
import os
import time
from pathlib import Path

import pytest

import htg
from htg.runner import WorkerPool, default_workers

TEMPLATE = Path(htg.__file__).parent / 'templates' / 'case_template.toml'


def fake_worker(worker_id, inbox, outbox, model_args, warm_up):
    """Speaks the worker protocol; 'solves' a point by returning a text.

    T = 999 kills the process (a crash); T = 998 sleeps, to test cancelling.
    """
    outbox.put(('state', worker_id, 'idle', ''))
    while True:
        message = inbox.get()
        if message is None:
            break
        task_id, case, chain = message
        for T, P in chain:
            if T == 999:
                os._exit(1)
            if T == 998:
                time.sleep(0.5)
            outbox.put(('point', worker_id, task_id, T, P, f'fake {worker_id}', 'test'))
        outbox.put(('task_done', worker_id, task_id))
    outbox.put(('state', worker_id, 'stopped', ''))


def case(**changes):
    changes.setdefault('aqueous_cutoff_C', None)
    return htg.load_case(TEMPLATE).replace(**changes)


def test_default_workers_positive():
    assert default_workers() >= 1


def test_all_points_of_all_runs_are_solved(tmp_path):
    with WorkerPool(3, worker_target=fake_worker) as pool:
        a = pool.submit(case(name='a', temperature_C=[300, 350, 400], pressure_bar=[240, 300]))
        b = pool.submit(case(name='b', temperature_C=[300, 400], pressure_bar=240, continuation=False),
                        output_dir=tmp_path / 'b')
        assert a.wait(20) and b.wait(20)
    assert a.status == b.status == 'done'
    assert set(a.outcomes) == set(a.case.grid()) and a.total == 6
    assert set(b.outcomes) == set(b.case.grid())
    assert (tmp_path / 'b' / 'summary.csv').exists() and (tmp_path / 'b' / 'case.toml').exists()
    assert [p[:2] for p in a.sweep_result().points] == a.case.grid()


def test_chain_runs_on_one_worker():
    with WorkerPool(2, worker_target=fake_worker) as pool:
        run = pool.submit(case(temperature_C=[300, 350, 400, 450], pressure_bar=240))
        assert run.wait(20)
    assert len({outcome for outcome, _ in run.outcomes.values()}) == 1


def test_cancel_drops_queued_pieces():
    with WorkerPool(1, worker_target=fake_worker) as pool:
        first = pool.submit(case(name='slow', temperature_C=[998, 998.5], pressure_bar=240,
                                 continuation=False))
        second = pool.submit(case(name='later', temperature_C=[300, 400], pressure_bar=240))
        pool.cancel(second)
        assert second.wait(5) and second.status == 'cancelled' and not second.outcomes
        assert first.wait(20) and first.status == 'done'


def test_crashed_worker_is_replaced():
    with WorkerPool(1, worker_target=fake_worker) as pool:
        crash = pool.submit(case(name='crash', temperature_C=[300, 999], pressure_bar=240))
        after = pool.submit(case(name='after', temperature_C=[300], pressure_bar=240))
        assert crash.wait(30) and after.wait(30)
    assert crash.outcomes[(999.0, 240.0)][0] == 'worker process stopped unexpectedly'
    assert after.status == 'done' and len(after.outcomes) == 1


@pytest.mark.slow
def test_real_workers_match_serial(base_model):
    """Two real workers (Julia) give the same results as solving in this process."""
    from conftest import FROZEN_DATA
    c = case(temperature_C=[300, 450], pressure_bar=[240, 300], continuation=False, aqueous_cutoff_C=390)
    with WorkerPool(2, {'data_dir': str(FROZEN_DATA)}, warm_up=False) as pool:
        run = pool.submit(c)
        assert run.wait(600)
    model = htg.Model.load(FROZEN_DATA)
    for (T, P), (result, _) in run.outcomes.items():
        assert not isinstance(result, str), result
        assert result.gibbs_J_per_mol_C == pytest.approx(model.solve(c, T, P).gibbs_J_per_mol_C, rel=1e-10)
