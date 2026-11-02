# Sweep 2 — best *constant* step size per instance, designed with the
# PEP_NumericalDesign optimizers instead of sampled on a grid.
#
# Same instance grid as `sweep_alpha.jl`, except that `samples_alpha` now plays
# the role of the multistart initial points: for every (T, ‖A‖, κ) and every
# start α₀ the constant policy is optimized with SOM (damped Newton on the
# frozen-certificate curvature) and/or FOM (normalized steepest descent). One
# row per (instance, start, optimizer) — aggregate offline with DataFrames:
#
#     using Arrow, DataFrames
#     df = DataFrame(Arrow.Table("Julia/results/<run>/merged.arrow"))
#     combine(groupby(df, [:T, :norm_A, :kappa]), :pep_best => minimum)
#
# Usage:
#   julia --project=Julia --threads=auto Julia/scripts/design_best_alpha.jl [--smoke] [--out DIR]

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

include(joinpath(@__DIR__, "..", "src", "PEPControl.jl"))
using .PEPControl
using Printf, Dates

const SMOKE = any(a -> a in ("--smoke", "-s"), ARGS)
const DEFAULT_ITERS = SMOKE ? 8 : 30
const DESIGN_ITERS = get(ENV, "PEP_CLUSTER_CONFIG", "0") == "1" ?
    parse(Int, get(ENV, "PEP_DESIGN_ITERS", string(DEFAULT_ITERS))) :
    DEFAULT_ITERS

const CFG = validate_common_cfg(cluster_overrides((
    T_lo = 1, T_hi = SMOKE ? 2 : 5,
    A_lo = 0.0, A_hi = 1.5, n_A = SMOKE ? 3 : 15,
    # multistart initial step sizes (kept small: each start costs many solves)
    alpha_lo = 0.05, alpha_hi = 1.5, n_alpha = SMOKE ? 2 : 4,
    # normalized strong-convexity grid κ = μ_f/L_f (linear, including κ = 0)
    kappa_lo = 0.0, kappa_hi = 1.0, n_kappa = SMOKE ? 2 : 8,
    optimizers = cluster_optimizers(SMOKE ? [:som] : [:som, :fom]),
    iters = DESIGN_ITERS,
    norm_B = 1.0, R_x = 1.0, ratio_radius = 1.0, L_f = 1.0,
    facial_reduction = :none,
    fr_rtol = 1e-9, fr_atol = 1e-11,
    cert_tol = 1e-7, gap_tol = 1e-6,
    fd_rel_step = 1e-4,
    batch_rows = 32,
)))

CFG.iters >= 1 || error("design iterations must be positive")

CFG.facial_reduction != :none && any(!=(:fom), CFG.optimizers) &&
    error("facial reduction supports only the FOM designer; set PEP_OPTIMIZERS=fom")

const RUN_DIR = let i = findfirst(==("--out"), ARGS)
    i === nothing ? joinpath(@__DIR__, "..", "results", "design_$(cfg_hash(CFG))") :
                    ARGS[i + 1]
end

const SAMPLES_T = cluster_horizons(CFG)
const SAMPLES_A = collect(range(CFG.A_lo, CFG.A_hi; length = CFG.n_A))
const SAMPLES_ALPHA = collect(range(CFG.alpha_lo, CFG.alpha_hi; length = CFG.n_alpha))
const SAMPLES_KAPPA = kappa_grid(CFG.kappa_lo, CFG.kappa_hi, CFG.n_kappa)

# The CECI `all` pipeline initializes its single FOM run from the best sampled
# alpha of the preceding sweep. Local runs retain the CFG multistart grid.
const USE_SWEEP_INITIALIZATION = get(ENV, "PEP_DESIGN_INIT_FROM_SWEEP", "0") == "1"
const FALLBACK_ALPHA = parse(Float64, get(ENV, "PEP_DESIGN_FALLBACK_ALPHA", "0.5"))
isfinite(FALLBACK_ALPHA) || error("PEP_DESIGN_FALLBACK_ALPHA must be finite")
const SWEEP_RESULTS = get(ENV, "PEP_SWEEP_RESULTS", "")
const SWEEP_INITIAL_ALPHAS = USE_SWEEP_INITIALIZATION ?
    best_sweep_alphas(SWEEP_RESULTS) :
    Dict{Tuple{Int, Float64, Float64}, Float64}()

function initial_alphas(T, normA, κ)
    USE_SWEEP_INITIALIZATION || return SAMPLES_ALPHA
    return [get(SWEEP_INITIAL_ALPHAS, (Int(T), Float64(normA), Float64(κ)),
                FALLBACK_ALPHA)]
end

const Row = @NamedTuple{T::Int, n_steps::Int, kappa::Float64, mu_f::Float64,
                        L_f::Float64, norm_A::Float64, norm_B::Float64,
                        R_x::Float64, R_u::Float64, alpha_init::Float64,
                        optimizer::String, alpha_best::Float64,
                        pep_best::Float64, n_solves::Int, status::String,
                        solve_time::Float64}

const KEY = (:T, :norm_A, :kappa, :alpha_init, :optimizer)

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
    nstarts = USE_SWEEP_INITIALIZATION ? 1 : length(SAMPLES_ALPHA)
    per_group = nstarts * length(CFG.optimizers)
    ntot = length(groups) * per_group
    @info "design_best_alpha" run_dir = RUN_DIR task_dir task_id = array.task_id task_count = array.task_count groups = length(groups) total_groups = length(all_groups) designs = ntot resumed = length(done) threads = Threads.nthreads() sweep_initialization = USE_SWEEP_INITIALIZATION sweep_groups = length(SWEEP_INITIAL_ALPHAS) fallback_alpha = FALLBACK_ALPHA

    finished = Threads.Atomic{Int}(0)
    t0 = time()

    pforeach(groups) do (T, normA, κ)
        starts = initial_alphas(T, normA, κ)
        todo = [(a, opt) for a in starts for opt in CFG.optimizers
                if (T, normA, κ, a, String(opt)) ∉ done]
        if isempty(todo)
            Threads.atomic_add!(finished, per_group)
            return
        end
        inst = ControlInstance(; T, norm_A = normA, norm_B = CFG.norm_B,
                               R_x = CFG.R_x, R_u = CFG.R_x * CFG.ratio_radius,
                               mu_f = κ * CFG.L_f, L_f = CFG.L_f)
        n = nsteps(inst)
        cp, cstatus = try_compile(inst)
        for (a0, opt) in todo
            t = time()
            r = if cp === nothing
                (alpha_best = NaN, pep_best = NaN, n_solves = 0,
                 status = cstatus)
            elseif opt == :fom && CFG.facial_reduction != :none
                # Optimize log(alpha), so every designed constant step
                # remains positive.
                best_positive_constant_alpha(
                    cp, n, a0; iters = CFG.iters,
                    facial_reduction = CFG.facial_reduction,
                    fr_rtol = CFG.fr_rtol, fr_atol = CFG.fr_atol,
                    cert_tol = CFG.cert_tol, gap_tol = CFG.gap_tol,
                    fd_rel_step = CFG.fd_rel_step)
            else
                best_constant_alpha(
                    cp, n, a0; optimizer = opt, iters = CFG.iters,
                    facial_reduction = CFG.facial_reduction,
                    fr_rtol = CFG.fr_rtol, fr_atol = CFG.fr_atol,
                    cert_tol = CFG.cert_tol, gap_tol = CFG.gap_tol,
                    fd_rel_step = CFG.fd_rel_step)
            end
            record!(writer, Row((inst.T, n, κ, inst.mu_f, inst.L_f, inst.norm_A,
                                 inst.norm_B, inst.R_x, inst.R_u, a0,
                                 String(opt), r.alpha_best, r.pep_best,
                                 r.n_solves, r.status, time() - t)))
        end
        k = Threads.atomic_add!(finished, per_group) + per_group
        @info @sprintf("[%5.1f%%] T=%d ‖A‖=%.3f κ=%.3g (%.1f s elapsed)",
                       100k / ntot, T, normA, κ, time() - t0)
    end

    nrows = close_writer!(writer)
    open(joinpath(task_dir, "_SUCCESS"), "w") do io
        println(io, "completed=", Dates.now())
    end
    path = array.task_count == 1 ? merge_parts(RUN_DIR) : nothing
    @info "design_best_alpha: done" task_id = array.task_id new_rows = nrows merged = path walltime = time() - t0
end

main()
