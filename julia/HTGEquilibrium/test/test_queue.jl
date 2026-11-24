# The worker pool (threads): scheduling with a stand-in chain solver (no
# solving), plus one real run (Python: tests/test_runner.py, which uses
# stand-in worker processes).
using Test
using HTGEquilibrium
using HTGEquilibrium: OrderedDict, start
isdefined(Main, :DataTestSetup) || include("data_setup.jl")
using Main.DataTestSetup: MODELS, REPO

const QUEUE_TEMPLATE = joinpath(REPO, "src", "htg", "templates", "case_template.toml")
queue_case(; changes...) = replace(load_case(QUEUE_TEMPLATE); aqueous_cutoff_C=nothing, changes...)

"""
Stand-in for solve_chain: "solves" a point by returning a text naming the task
that solved it. T = 999 throws (an internal error); T = 998 sleeps, to test
cancelling and the status display.
"""
function fake_chain(model, case, chain, on_point; stop=() -> false)
    for (T, P) in chain
        stop() && return
        T == 999 && error("crash in the equation of state")
        T == 998 && sleep(0.5)
        on_point(T, P, "fake $(objectid(current_task()))", "test")
    end
end

fake_pool(f, n; kwargs...) = WorkerPool(f, n; model=MODELS["base"], warm_up=false, chain_solver=fake_chain, kwargs...)

@testset "queue" begin
    @test default_workers() >= 1

    @testset "all points of all runs are solved" begin
        dir = mktempdir()
        a, b = fake_pool(3) do pool
            a = submit(pool, queue_case(name="a", temperature_C=[300, 350, 400], pressure_bar=[240, 300]))
            b = submit(pool, queue_case(name="b", temperature_C=[300, 400], pressure_bar=240, continuation=false),
                       joinpath(dir, "b"))
            @test wait(a, 20) && wait(b, 20)
            a, b
        end
        @test a.status == b.status == "done"
        @test Set(keys(a.outcomes)) == Set(grid(a.case)) && a.total == 6 && a.done_points == 6
        @test a.failed_points == 6   # the stand-in returns texts, which count as failures
        @test Set(keys(b.outcomes)) == Set(grid(b.case))
        @test isfile(joinpath(dir, "b", "summary.csv")) && isfile(joinpath(dir, "b", "case.toml"))
        @test last(b.log) == "results saved to $(joinpath(dir, "b"))"
        @test [p[1:2] for p in sweep_result(a).points] == grid(a.case)
        @test a.name == "a" && a.finished_running && a.started !== nothing && a.finished >= a.started
    end

    @testset "a chain runs on one worker" begin
        run = fake_pool(2) do pool
            run = submit(pool, queue_case(temperature_C=[300, 350, 400, 450], pressure_bar=240))
            @test wait(run, 20)
            run
        end
        @test length(Set(first.(values(run.outcomes)))) == 1
    end

    @testset "cancel drops queued pieces" begin
        fake_pool(1) do pool
            first_run = submit(pool, queue_case(name="slow", temperature_C=[998, 998.5], pressure_bar=240,
                                                continuation=false))
            second = submit(pool, queue_case(name="later", temperature_C=[300, 400], pressure_bar=240))
            cancel(pool, second)
            @test wait(second, 5) && second.status == "cancelled" && isempty(second.outcomes)
            @test last(second.log) == "cancelled"
            @test wait(first_run, 20) && first_run.status == "done"
        end
    end

    @testset "an internal error fails the rest of the chain; the worker goes on" begin
        crash, after = fake_pool(1) do pool
            crash = submit(pool, queue_case(name="crash", temperature_C=[300, 999], pressure_bar=240))
            after = submit(pool, queue_case(name="after", temperature_C=[300], pressure_bar=240))
            @test wait(crash, 30) && wait(after, 30)
            crash, after
        end
        @test startswith(crash.outcomes[(999.0, 240.0)][1], "internal error: crash in the equation of state")
        @test crash.status == "done" && crash.failed_points == 2
        @test after.status == "done" && length(after.outcomes) == 1
    end

    @testset "status, set_workers, callbacks, restart, shutdown" begin
        points, done = Any[], Any[]
        pool = start(WorkerPool(2; model=MODELS["base"], warm_up=false, chain_solver=fake_chain,
                                on_point=(run, T, P, outcome, start) -> push!(points, (run.name, T, P, start)),
                                on_done=run -> push!(done, run.name)))
        try
            slow = submit(pool, queue_case(name="slow", temperature_C=[998, 998.5, 998.7], pressure_bar=240))
            sleep(0.2)
            s = status(pool)
            @test s["target"] == 2 && s["queued_pieces"] == 0 && length(s["workers"]) == 2
            busy = [w for w in s["workers"] if w["state"] == "busy"]
            @test length(busy) == 1 && busy[1]["case"] == "slow" && busy[1]["point"] == (998.0, 240.0)
            # a new model is used for runs submitted from now on; the running one keeps its model
            @test restart(pool, MODELS["acetone"]) == false
            later = submit(pool, queue_case(name="later", temperature_C=[300], pressure_bar=240))
            @test later.model === MODELS["acetone"] && slow.model === MODELS["base"]
            @test wait(slow, 20) && wait(later, 20)
            @test restart(pool, MODELS["base"]) == true
            @test Set(done) == Set(["slow", "later"]) && length(points) == 4
            @test (("slow", 998.0, 240.0, "test") in points)
            set_workers(pool, 4)
            @test status(pool)["target"] == 4 && count(w -> w.state != "stopped", values(pool.workers)) == 4
            set_workers(pool, 1)
            sleep(0.2)   # idle extra workers stop
            @test length(status(pool)["workers"]) == 1
            many = [submit(pool, queue_case(name="m$k", temperature_C=[300, 350], pressure_bar=240)) for k in 1:3]
            @test all(r -> wait(r, 20), many)
        finally
            shutdown(pool)
        end
        @test all(w -> w.state == "stopped", values(pool.workers))
    end

    @testset "warm-up: one worker solves the template first" begin
        pool = start(WorkerPool(2; model=MODELS["base"], chain_solver=fake_chain))   # warm_up=true
        try
            states = [w["state"] for w in status(pool)["workers"]]
            @test count(==("warming"), states) <= 1
            run = submit(pool, queue_case(temperature_C=[300, 350], pressure_bar=240, continuation=false))
            @test wait(run, 600) && run.status == "done" && run.done_points == 2
        finally
            shutdown(pool)
        end
        @test count(w -> w.warm_up, values(pool.workers)) == 1
    end

    @testset "real solves on $(min(2, Threads.nthreads())) worker threads equal serial solves" begin
        model = MODELS["base"]
        c = queue_case(temperature_C=[300, 450], pressure_bar=[240, 300], continuation=false, aqueous_cutoff_C=390)
        run = WorkerPool(2; model, warm_up=false) do pool
            run = submit(pool, c)
            @test wait(run, 600)
            run
        end
        for ((T, P), (result, start)) in run.outcomes
            @test result isa Result
            @test start == "feasible LP"
            @test result.gibbs_J_per_mol_C == solve(model, c, T, P).gibbs_J_per_mol_C
        end
        @test summary(sweep_result(run))["status"] == fill("solved", length(grid(c))) == fill("solved", 6)
        tables = live_tables(run)
        @test tables.label == "my_case (this session, #$(run.id))" && nrows(tables.summary) == 6
    end
end
