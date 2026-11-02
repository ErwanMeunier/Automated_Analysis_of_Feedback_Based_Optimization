# Sweep 1 — worst-case PEP value on the grid
#     samples_A × samples_T × samples_alpha × samples_kappa
#
# Instances are grouped by (T, ‖A‖, κ) so the PEP is compiled once per group
# and reused for every step size; groups are distributed dynamically over the
# Julia threads. Every batch of rows is flushed to a self-contained Arrow shard
# under `results/<run>/`, so the run is resumable and readable offline:
#
#     using Arrow, DataFrames
#     df = DataFrame(Arrow.Table("Julia/results/<run>/merged.arrow"))
#
# Usage:
#   julia --project=Julia --threads=auto Julia/scripts/sweep_alpha.jl [--smoke] [--out DIR]

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

include(joinpath(@__DIR__, "..", "src", "PEPControl.jl"))
using .PEPControl
using Printf, Dates

const SMOKE = any(a -> a in ("--smoke", "-s"), ARGS)

const CFG = validate_common_cfg(cluster_overrides((
    # horizons T (the method performs T+1 updates)
    T_lo = 2, T_hi = SMOKE ? 10 : 20,
    # ‖A‖ grid (linear)
    A_lo = 0.0, A_hi = 1.0, n_A = SMOKE ? 3 : 10,
    # constant step-size grid (linear)
    alpha_lo = 1.0e-3, alpha_hi = 2.0, n_alpha = SMOKE ? 4 : 100,
    # normalized strong-convexity grid κ = μ_f/L_f (linear, including κ = 0)
    kappa_lo = 0.0, kappa_hi = 1.0, n_kappa = SMOKE ? 2 : 10,
    # fixed dimensionless constants
    norm_B = 1.0, R_x = 1.0, ratio_radius = 1.0, L_f = 1.0,
    facial_reduction = :none,
    fr_rtol = 1e-9, fr_atol = 1e-11,
    cert_tol = 1e-7, gap_tol = 1e-6,
    batch_rows = 64,
)))

const RUN_DIR = let i = findfirst(==("--out"), ARGS)
    i === nothing ? joinpath(@__DIR__, "..", "results", "sweep_$(cfg_hash(CFG))") :
                    ARGS[i + 1]
end

const SAMPLES_T = cluster_horizons(CFG)
const SAMPLES_A = collect(range(CFG.A_lo, CFG.A_hi; length = CFG.n_A))
const SAMPLES_ALPHA = collect(range(CFG.alpha_lo, CFG.alpha_hi; length = CFG.n_alpha))
const SAMPLES_KAPPA = kappa_grid(CFG.kappa_lo, CFG.kappa_hi, CFG.n_kappa)

const Row = @NamedTuple{T::Int, n_steps::Int, kappa::Float64, mu_f::Float64,
                        L_f::Float64, norm_A::Float64, norm_B::Float64,
                        R_x::Float64, R_u::Float64, alpha::Float64,
                        pep::Float64, status::String, solve_time::Float64}

const KEY = (:T, :norm_A, :kappa, :alpha)

function main()
    mkpath(RUN_DIR)
    array = array_context()
    task_dir = array_task_dir(RUN_DIR, array.task_id, array.task_count)
    mkpath(task_dir)
    array.task_id == 0 && write_config(
        RUN_DIR, merge(CFG, (T_values = SAMPLES_T,
                             kappa_definition = "mu_f/L_f")))
    done = completed_keys(task_dir, KEY)
    writer = ResultWriter(Row, task_dir; batch = CFG.batch_rows)

    all_groups = [(T, normA, κ) for T in SAMPLES_T, normA in SAMPLES_A,
                  κ in SAMPLES_KAPPA] |> vec
    groups = array_partition(all_groups, array.task_id, array.task_count)
    ntot = length(groups) * length(SAMPLES_ALPHA)
    @info "sweep_alpha" run_dir = RUN_DIR task_dir task_id = array.task_id task_count = array.task_count groups = length(groups) total_groups = length(all_groups) instances = ntot resumed = length(done) threads = Threads.nthreads()

    finished = Threads.Atomic{Int}(0)
    t0 = time()

    pforeach(groups) do (T, normA, κ)
        todo = [a for a in SAMPLES_ALPHA if (T, normA, κ, a) ∉ done]
        if isempty(todo)
            Threads.atomic_add!(finished, length(SAMPLES_ALPHA))
            return
        end
        inst = ControlInstance(; T, norm_A = normA, norm_B = CFG.norm_B,
                               R_x = CFG.R_x, R_u = CFG.R_x * CFG.ratio_radius,
                               mu_f = κ * CFG.L_f, L_f = CFG.L_f)
        n = nsteps(inst)
        cp, cstatus = try_compile(inst)
        for a in todo
            t = time()
            W, status = cp === nothing ? (NaN, cstatus) :
                        pep_at_constant_alpha(
                            cp, n, a;
                            facial_reduction = CFG.facial_reduction,
                            fr_rtol = CFG.fr_rtol, fr_atol = CFG.fr_atol,
                            cert_tol = CFG.cert_tol, gap_tol = CFG.gap_tol)
            record!(writer, Row((inst.T, n, κ, inst.mu_f, inst.L_f, inst.norm_A,
                                 inst.norm_B, inst.R_x, inst.R_u, a, W, status,
                                 time() - t)))
        end
        k = Threads.atomic_add!(finished, length(SAMPLES_ALPHA)) + length(SAMPLES_ALPHA)
        @info @sprintf("[%5.1f%%] T=%d ‖A‖=%.3f κ=%.3g (%.1f s elapsed)",
                       100k / ntot, T, normA, κ, time() - t0)
    end

    nrows = close_writer!(writer)
    open(joinpath(task_dir, "_SUCCESS"), "w") do io
        println(io, "completed=", Dates.now())
    end
    path = array.task_count == 1 ? merge_parts(RUN_DIR) : nothing
    @info "sweep_alpha: done" task_id = array.task_id new_rows = nrows merged = path walltime = time() - t0
end

main()
