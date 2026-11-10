"""Running cases on a pool of worker processes.

Each worker process loads the model once and stays alive, so the Julia start-up
(about a minute for the first solve) is paid once per worker, not per run.

Every submitted case is split into independent pieces (see Case.chains): single
points, or continuation chains that must run in order. All pieces of all runs go
into one first-in, first-out queue; an idle worker takes the next piece. One large
case spreads over all workers, and with several cases the one submitted first
finishes first.

    pool = WorkerPool(workers=4)
    run = pool.submit(case, output_dir='results/my_case')
    run.wait()
    run.sweep_result().summary()
    pool.shutdown()
"""
import itertools
import multiprocessing as mp
import os
import queue
import threading
import time
import traceback
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path

from .equilibrium import SolveError, apply_phase_rules
from .results import SweepResult

GB = 1024 ** 3
WORKER_MEMORY_GB = 1.5  # measured peak about 1.4 GB per warm worker


# ---------------------------------------------------------------- solving a chain
def solve_chain(model, case, chain, on_point):
    """Solve the (T, P) points of ``chain`` in order; call on_point(T, P, outcome, start).

    outcome is a Result or an error message. With continuation each point starts
    from the previous successful solution when the selected (phase, species)
    pairs are unchanged; otherwise from a fresh feasible LP.
    """
    previous = None
    for T, P in chain:
        point = case.at(T, P)
        start = 'feasible LP'
        if case.continuation and previous is not None:
            phases = apply_phase_rules(point, model.registry)[0]
            signature = {(p, s) for p, names in phases.items() for s in names}
            if set(previous.amounts) == signature:
                point.initial_amounts = dict(previous.amounts)
                start = f'solution at {previous.temperature_C:g} °C, {previous.pressure_bar:g} bar'
            else:
                point.initial_amounts = None
                start = 'feasible LP (phase selection changed)'
        try:
            outcome = model.solve(point)
            previous = outcome
        except SolveError as exc:
            outcome = str(exc)
        on_point(T, P, outcome, start)


def default_workers():
    """Number of worker processes: min(CPU cores - 1, available memory / 1.5 GB), at least 1."""
    cores = max(1, (os.cpu_count() or 2) - 1)
    try:
        import psutil
        memory = int(psutil.virtual_memory().available / (WORKER_MEMORY_GB * GB))
    except Exception:
        memory = cores
    return max(1, min(cores, memory))


# ---------------------------------------------------------------- runs
@dataclass
class Run:
    """A submitted case. Updated by the pool while it runs."""
    id: int
    case: object
    output_dir: Path | None
    total: int
    status: str = 'queued'  # queued, running, done, cancelled
    outcomes: dict = field(default_factory=dict)  # (T, P) -> (Result | str, start)
    log: list = field(default_factory=list)
    submitted: float = field(default_factory=time.time)
    started: float | None = None
    finished: float | None = None
    tasks_left: int = 0
    _event: threading.Event = field(default_factory=threading.Event, repr=False)

    @property
    def name(self):
        return self.case.name

    @property
    def done_points(self):
        return len(self.outcomes)

    @property
    def failed_points(self):
        return sum(1 for outcome, _ in self.outcomes.values() if isinstance(outcome, str))

    @property
    def finished_running(self):
        return self.status in ('done', 'cancelled')

    def wait(self, timeout=None):
        return self._event.wait(timeout)

    def sweep_result(self):
        """Results so far, in grid order (points not solved yet are left out)."""
        points = [(T, P, self.outcomes[(T, P)][0]) for T, P in self.case.grid()
                  if (T, P) in self.outcomes]
        return SweepResult(self.case, points)


@dataclass
class _Task:
    id: int
    run: Run
    chain: list


@dataclass
class _Worker:
    id: int
    process: object
    inbox: object
    state: str = 'starting'  # starting, warming, idle, busy, stopping, stopped, crashed
    task: _Task | None = None
    retire: bool = False
    message: str = ''


# ---------------------------------------------------------------- worker process
def _worker_main(worker_id, inbox, outbox, model_args, warm_up):
    try:
        from .case import load_case
        from .model import Model
        outbox.put(('state', worker_id, 'starting', 'loading databases'))
        model = Model.load(**model_args)
        if warm_up:
            outbox.put(('state', worker_id, 'warming', 'compiling the equation of state'))
            try:
                template = Path(__file__).parent / 'templates' / 'case_template.toml'
                model.solve(load_case(template).at(300.0, 240.0))
            except Exception:
                pass
    except Exception:
        outbox.put(('state', worker_id, 'crashed', traceback.format_exc()[-500:]))
        return
    outbox.put(('state', worker_id, 'idle', ''))
    while True:
        message = inbox.get()
        if message is None:
            break
        task_id, case, chain = message

        def on_point(T, P, outcome, start):
            outbox.put(('point', worker_id, task_id, T, P, outcome, start))
        try:
            solve_chain(model, case, chain, on_point)
        except Exception:
            outbox.put(('task_error', worker_id, task_id, traceback.format_exc()[-800:]))
        outbox.put(('task_done', worker_id, task_id))
    outbox.put(('state', worker_id, 'stopped', ''))


# ---------------------------------------------------------------- pool
class WorkerPool:
    """Persistent worker processes with one FIFO queue of chains.

    workers     number of processes (default: default_workers())
    model_args  arguments of Model.load in the workers (data_dir, species_file, ...)
    warm_up     solve a small case right after start so the first real point is fast
    on_point    optional callback(run, T, P, outcome, start), called from a background thread
    on_done     optional callback(run)
    worker_target  process function (tests replace it with a stand-in without Julia)
    """

    def __init__(self, workers=None, model_args=None, warm_up=True, on_point=None, on_done=None,
                 worker_target=None):
        self.worker_target = worker_target or _worker_main
        self.target = int(workers or default_workers())
        self.model_args = dict(model_args or {})
        self.warm_up = warm_up
        self.on_point, self.on_done = on_point, on_done
        self._ctx = mp.get_context('spawn')  # Julia must not be forked
        self._outbox = self._ctx.Queue()
        self._workers = {}
        self._pending = deque()
        self._tasks = {}
        self.runs = []
        self._ids = itertools.count(1)
        self._lock = threading.RLock()
        self._stop = threading.Event()
        self._thread = None

    # ------------------------------------------------------------ public API
    def start(self):
        if self._thread is None:
            with self._lock:
                for _ in range(self.target):
                    self._spawn()
            self._thread = threading.Thread(target=self._loop, name='htg-dispatcher', daemon=True)
            self._thread.start()
        return self

    def __enter__(self):
        return self.start()

    def __exit__(self, *exc):
        self.shutdown()

    def submit(self, case, output_dir=None):
        """Queue ``case``; returns a Run. Results are saved to output_dir when it finishes."""
        case.validate()
        chains = case.chains()
        with self._lock:
            run = Run(id=next(self._ids), case=case,
                      output_dir=Path(output_dir) if output_dir else None,
                      total=sum(len(c) for c in chains), tasks_left=len(chains))
            for chain in chains:
                task = _Task(id=next(self._ids), run=run, chain=chain)
                self._tasks[task.id] = task
                self._pending.append(task)
            self.runs.append(run)
        self.start()
        return run

    def cancel(self, run):
        """Drop the run's queued pieces; pieces already running finish their current chain."""
        with self._lock:
            if run.finished_running:
                return
            kept = deque()
            for task in self._pending:
                if task.run is run:
                    run.tasks_left -= 1
                    self._tasks.pop(task.id, None)
                else:
                    kept.append(task)
            self._pending = kept
            run.log.append('cancelled')
            run.status = 'cancelling'
            if run.tasks_left == 0:
                self._finish(run, 'cancelled')

    def set_workers(self, n):
        """Change the number of worker processes (idle extra workers stop)."""
        with self._lock:
            self.target = max(1, int(n))
            active = [w for w in self._workers.values() if w.state not in ('stopped', 'crashed') and not w.retire]
            for _ in range(self.target - len(active)):
                self._spawn()
            for worker in active[self.target:]:
                worker.retire = True

    def status(self):
        """Worker states and queue length (for display)."""
        with self._lock:
            workers = [{'id': w.id, 'state': w.state, 'message': w.message,
                        'case': w.task.run.name if w.task else ''}
                       for w in self._workers.values() if w.state not in ('stopped',)]
            return {'workers': workers, 'queued_pieces': len(self._pending), 'target': self.target}

    def shutdown(self, wait=True):
        self._stop.set()
        with self._lock:
            workers = list(self._workers.values())
        for worker in workers:
            if worker.process.is_alive():
                try:
                    worker.inbox.put(None)
                except Exception:
                    pass
        if wait:
            for worker in workers:
                worker.process.join(timeout=5)
        for worker in workers:
            if worker.process.is_alive():
                worker.process.terminate()
        if self._thread is not None:
            self._thread.join(timeout=2)

    # ------------------------------------------------------------ internals
    def _spawn(self):
        worker_id = next(self._ids)
        inbox = self._ctx.Queue()
        process = self._ctx.Process(target=self.worker_target, name=f'htg-worker-{worker_id}', daemon=True,
                                    args=(worker_id, inbox, self._outbox, self.model_args, self.warm_up))
        process.start()
        self._workers[worker_id] = _Worker(worker_id, process, inbox)

    def _loop(self):
        while not self._stop.is_set():
            try:
                message = self._outbox.get(timeout=0.2)
            except queue.Empty:
                message = None
            except (EOFError, OSError):
                break
            with self._lock:
                if message is not None:
                    self._handle(message)
                self._check_workers()
                self._dispatch()

    def _handle(self, message):
        kind, worker_id = message[0], message[1]
        worker = self._workers.get(worker_id)
        if worker is None:
            return
        if kind == 'state':
            _, _, state, text = message
            if worker.state != 'busy' or state in ('crashed', 'stopped'):
                worker.state, worker.message = state, text
        elif kind == 'point':
            _, _, task_id, T, P, outcome, start = message
            task = self._tasks.get(task_id)
            if task is None:
                return
            run = task.run
            run.outcomes[(T, P)] = (outcome, start)
            text = 'solved' if not isinstance(outcome, str) else f'FAILED: {outcome}'
            run.log.append(f'{T:g} °C, {P:g} bar: {text} (start: {start})')
            if self.on_point:
                self._callback(self.on_point, run, T, P, outcome, start)
        elif kind == 'task_error':
            task = self._tasks.get(message[2])
            if task is not None:
                self._fail_rest(task, f'internal error: {message[3]}')
        elif kind == 'task_done':
            task = self._tasks.pop(message[2], None)
            worker.task, worker.state = None, 'idle'
            if task is not None:
                self._task_finished(task)

    def _check_workers(self):
        for worker in list(self._workers.values()):
            if worker.state in ('stopped', 'crashed') or worker.process.is_alive():
                continue
            # The process died (e.g. a crash inside the equation of state).
            worker.state, worker.message = 'crashed', f'exit code {worker.process.exitcode}'
            task, worker.task = worker.task, None
            if task is not None:
                self._tasks.pop(task.id, None)
                self._fail_rest(task, 'worker process stopped unexpectedly')
                self._task_finished(task)
            if not worker.retire and not self._stop.is_set():
                self._spawn()

    def _dispatch(self):
        for worker in self._workers.values():
            if worker.state != 'idle':
                continue
            if worker.retire:
                worker.state = 'stopping'
                worker.inbox.put(None)
                continue
            if not self._pending:
                break
            task = self._pending.popleft()
            run = task.run
            if run.status == 'queued':
                run.status, run.started = 'running', time.time()
            worker.task, worker.state = task, 'busy'
            worker.inbox.put((task.id, run.case, task.chain))

    def _fail_rest(self, task, reason):
        for T, P in task.chain:
            if (T, P) not in task.run.outcomes:
                task.run.outcomes[(T, P)] = (reason, '')
                task.run.log.append(f'{T:g} °C, {P:g} bar: FAILED: {reason}')

    def _task_finished(self, task):
        run = task.run
        run.tasks_left -= 1
        if run.tasks_left <= 0 and not run.finished_running:
            self._finish(run, 'cancelled' if run.status == 'cancelling' else 'done')

    def _finish(self, run, status):
        run.status, run.finished = status, time.time()
        if run.output_dir is not None and run.outcomes:
            try:
                run.sweep_result().save(run.output_dir)
                run.log.append(f'results saved to {run.output_dir}')
            except Exception as exc:
                run.log.append(f'saving results failed: {exc}')
        if self.on_done:
            self._callback(self.on_done, run)
        run._event.set()

    def _callback(self, function, *args):
        try:
            function(*args)
        except Exception:
            traceback.print_exc()
