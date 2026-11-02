# PEPControl — PEP models for the LTI control design problem, grid definitions,
# a work-stealing thread pool, and crash-safe Arrow result storage.
#
# The Arrow shards written here are read back with
#   using Arrow, DataFrames; DataFrame(Arrow.Table("results/<run>/merged.arrow"))
# or with `PEPControl.load_results("results/<run>")`.

module PEPControl

include(joinpath(@__DIR__, "..", "PEP_NumericalDesign", "src", "PEPDesign.jl"))
using .PEPDesign

using Arrow, DataFrames, Dates, Printf

export ControlInstance, nsteps, kappa, build_pep_control, compile_control,
       try_compile, pep_at_constant_alpha, best_constant_alpha,
       best_positive_constant_alpha, analytic_constant_alpha,
       ResultWriter, record!, flush_writer!, close_writer!,
       load_results, load_array_results, merge_parts, merge_array_parts,
       completed_keys, write_config, array_context, array_partition,
       array_task_dir, cluster_overrides, cluster_optimizers,
       cluster_horizons, validate_common_cfg, best_sweep_alphas,
       pforeach, cfg_hash, log_grid, kappa_grid

# ── problem instance ─────────────────────────────────────────────────────────

"""
    ControlInstance(; T, norm_A, norm_B, R_x, R_u, mu_f, L_f)

One PEP-control instance: horizon `T` (the method performs `T+1` updates),
operator norms `‖A‖ ≤ norm_A`, `‖B‖ ≤ norm_B`, initial radii `R_x`, `R_u`,
and the objective class `F_{mu_f, L_f}`.
"""
Base.@kwdef struct ControlInstance
    T::Int
    norm_A::Float64
    norm_B::Float64 = 1.0
    R_x::Float64 = 1.0
    R_u::Float64 = 1.0
    mu_f::Float64
    L_f::Float64 = 1.0
end

"""Number of step-size coefficients of the compiled PEP (= horizon + 1)."""
nsteps(inst::ControlInstance) = inst.T + 1

"""Normalized strong-convexity ratio `mu_f / L_f` in `[0, 1]`."""
kappa(inst::ControlInstance) = inst.mu_f / inst.L_f

function _validate_instance(inst::ControlInstance)
    inst.T >= 0 || error("T must be nonnegative")
    isfinite(inst.L_f) && inst.L_f > 0 ||
        error("L_f must be finite and positive")
    isfinite(inst.mu_f) && 0 <= inst.mu_f <= inst.L_f ||
        error("need 0 <= mu_f <= L_f (equivalently 0 <= kappa = mu_f/L_f <= 1)")
    return inst
end

"""Create the exact smooth convex class associated with `mu_f / L_f`."""
function _objective_class(m::PEPModel, inst::ControlInstance)
    if iszero(inst.mu_f)
        return SmoothConvexFunction(m; L = inst.L_f)
    elseif inst.mu_f == inst.L_f
        # The generic F_{mu,L} interpolation divides by L-mu. At equality,
        # the class consists exactly of quadratics with Hessian L*I.
        return SmoothStronglyConvexQuadraticFunction(
            m; μ = inst.mu_f, L = inst.L_f)
    end
    return SmoothStronglyConvexFunction(m; μ = inst.mu_f, L = inst.L_f)
end

# ── PEP model ────────────────────────────────────────────────────────────────

"""
    build_pep_control(inst) -> (model, alpha)

PEP model of the control-design loop: LTI dynamics `x⁺ = A x + B u`, inverse
gradient mapping `g_inv = A' g_inv + g`, and control update
`u⁺ = u − α_t B' g_inv_t`. Maximizes `f(x_T) − f_*`.
"""

function build_pep_control(inst::ControlInstance)
    _validate_instance(inst)
    T = nsteps(inst)
    m = PEPModel()
    alpha = coeffs!(m, T)
    f = _objective_class(m, inst)
    x_star = point!(m; name = "x_star")
    g_star, f_star = oracle!(f, x_star)
    g_star_inv = point!(m; name = "g_star_inv")

    u_star = point!(m; name = "u_star")
    A = LinearOperator(m; L = inst.norm_A, name = "A")
    B = LinearOperator(m; L = inst.norm_B, name = "B")

    add_eq!(m, sqnorm(g_star_inv - adjoint_oracle!(A, g_star_inv)[1] - g_star))  # enforce stationarity at x_star
    add_eq!(m, sqnorm(adjoint_oracle!(B, g_star_inv)[1]))  # enforce steady state at x_star

    g_inv = Vector{PointExpr}(undef, T)
    x = Vector{PointExpr}(undef, T + 1)
    u = Vector{PointExpr}(undef, T + 1)
    g = Vector{PointExpr}(undef, T)
    f_vec = Vector{Any}(undef, T)

    x[1] = point!(m; name = "x1")
    u[1] = point!(m; name = "u1")
    add_le!(m, sqnorm(x[1] - x_star) - inst.R_x^2)
    add_le!(m, sqnorm(u[1] - u_star) - inst.R_u^2)

    # steady state: x* = A x* + B u*
    add_eq!(m, sqnorm(x_star - gradient!(A, x_star) - gradient!(B, u_star)))

    for t in 1:T
        g_inv[t] = point!(m; name = "g_inv_$t")
        x[t + 1] = gradient!(A, x[t]) + gradient!(B, u[t])
        u[t + 1] = u[t] - alpha[t] * adjoint_oracle!(B, g_inv[t])[1]
        g[t], f_vec[t] = oracle!(f, x[t])
        add_eq!(m, sqnorm(g_inv[t] - adjoint_oracle!(A, g_inv[t])[1] - g[t]))
    end
    u_TM = point!(m; name = "u_TM")
    add_eq!(m, sqnorm(u_TM -gradient!(A, u_TM) - gradient!(B, u[T])))
    _ , fGUT = oracle!(f, u_TM) 
    objective_max!(m, fGUT - f_star)
    return m, alpha
end

"""Compile the instance once; the result is reusable for every step size."""
compile_control(inst::ControlInstance) = compile(build_pep_control(inst)[1])

"""
    try_compile(inst) -> (cp, status)

Non-throwing `compile_control`: a degenerate or unsupported instance yields
`(nothing, "ERROR: …")` so that a sweep records it instead of dying.
"""
function try_compile(inst::ControlInstance)
    try
        return compile_control(inst), "OK"
    catch err
        return nothing, "ERROR: " * first(sprint(showerror, err), 120)
    end
end

"""
    pep_at_constant_alpha(cp, n, alpha; facial_reduction=:none) -> (W, status)

Worst case of the compiled PEP under the constant policy `α_t ≡ alpha`.
Returns `(NaN, "ERROR: …")` instead of throwing when the solve fails.
With `facial_reduction=:explicit`, an uncertified value is returned as `NaN`
rather than being mistaken for a finite worst case.
"""
function pep_at_constant_alpha(cp::CompiledPEP, n::Int, alpha::Real;
                               facial_reduction::Symbol = :none,
                               fr_rtol::Float64 = 1e-9,
                               fr_atol::Float64 = 1e-11,
                               cert_tol::Float64 = 1e-7,
                               gap_tol::Float64 = 1e-6)
    try
        sol = solve_pep(cp, fill(Float64(alpha), n); warn = false,
                        facial_reduction, fr_rtol, fr_atol, cert_tol, gap_tol)
        if facial_reduction != :none && !sol.diagnostics.certified
            diagnostics = sol.diagnostics
            return NaN, @sprintf(
                "UNCERTIFIED: %s gap=%.3e scalar=%.3e cone=%.3e",
                string(sol.status), diagnostics.relative_gap,
                diagnostics.max_scalar_violation,
                diagnostics.max_cone_violation)
        end
        return sol.obj, string(sol.status)
    catch err
        return NaN, "ERROR: " * first(sprint(showerror, err), 120)
    end
end

"""
    best_constant_alpha(cp, n, alpha0; optimizer, iters) -> NamedTuple

Design the best constant step size from the start `alpha0` with the
PEP_NumericalDesign optimizers (`:som` = damped Newton on the
frozen-certificate curvature, `:fom` = normalized steepest descent).
"""
function best_constant_alpha(cp::CompiledPEP, n::Int, alpha0::Real;
                             optimizer::Symbol = :som, iters::Int = 30,
                             facial_reduction::Symbol = :none,
                             fr_rtol::Float64 = 1e-9,
                             fr_atol::Float64 = 1e-11,
                             cert_tol::Float64 = 1e-7,
                             gap_tol::Float64 = 1e-6,
                             fd_rel_step::Float64 = 1e-4)
    try
        facial_reduction != :none && optimizer != :fom &&
            error("facial reduction requires optimizer=:fom")
        dp = DesignProblem(cp, ConstantPolicy(), n;
                           facial_reduction, fr_rtol, fr_atol, cert_tol,
                           gap_tol, fd_rel_step)
        tr = optimizer === :som ? design_som(dp, [Float64(alpha0)]; iters) :
             optimizer === :fom ? design_fom(dp, [Float64(alpha0)]; iters,
                                             steps = t -> 0.1 / sqrt(t)) :
             error("unknown optimizer :$optimizer (use :som or :fom)")
        ω, W = best_point(tr)
        isfinite(ω[1]) && isfinite(W) ||
            error("designer produced no finite certified point")
        return (alpha_best = ω[1], pep_best = W, n_solves = tr.nsolves,
                status = "OK")
    catch err
        return (alpha_best = NaN, pep_best = NaN, n_solves = 0,
                status = "ERROR: " * first(sprint(showerror, err), 120))
    end
end

"""
    analytic_constant_alpha(norm_A, L_f, norm_B=1) -> Float64

Analytic constant-step bound used as the design initializer. The usual
`norm_B = 1` formula scales as `1 / norm_B^2`, because the plant-to-state map
is linear in `B`. Returns `NaN` outside its stable/identifiable domain.
"""
function analytic_constant_alpha(norm_A::Real, L_f::Real, norm_B::Real = 1.0)
    a, L, b = Float64(norm_A), Float64(L_f), Float64(norm_B)
    (!isfinite(a) || !isfinite(L) || !isfinite(b) ||
     a >= 1 || L <= 0 || b <= 0) && return NaN
    P = 1 / (1 - a^2)
    G = 1 / (1 - a)
    L_G = L * G^2
    L_tilde = L * G
    zeta = 2 * abs(a * P * G)
    omega = P * G^2
    denominator = L_G + L_tilde * (zeta + sqrt(zeta^2 + 4omega))
    return denominator > 0 ? 2 / (b^2 * denominator) : NaN
end

"""
    best_positive_constant_alpha(cp, n, alpha0; iters=30, ...) -> NamedTuple

Optimize the positive constant policy `alpha_t = exp(omega)` with FOM. This
is the common implementation used by the laptop and CECI design experiments.
With explicit facial reduction, every finite value is certified in the
original unreduced PEP before it can enter the design trace.
"""
function best_positive_constant_alpha(cp::CompiledPEP, n::Int, alpha0::Real;
                                      iters::Int = 30,
                                      facial_reduction::Symbol = :none,
                                      fr_rtol::Float64 = 1e-9,
                                      fr_atol::Float64 = 1e-11,
                                      cert_tol::Float64 = 1e-7,
                                      gap_tol::Float64 = 1e-6,
                                      fd_rel_step::Float64 = 1e-4)
    try
        isfinite(alpha0) && alpha0 > 0 ||
            error("alpha0 must be finite and positive")
        policy = FunctionPolicy(
            (omega, _, _) -> exp(omega[1]), 1; name = "positive_constant")
        problem = DesignProblem(cp, policy, n;
                                facial_reduction, fr_rtol, fr_atol, cert_tol,
                                gap_tol, fd_rel_step)
        trace = design_fom(problem, [log(Float64(alpha0))]; iters,
                           steps = iteration -> 0.1 / sqrt(iteration))
        omega, W = best_point(trace)
        alpha = exp(omega[1])
        isfinite(alpha) && isfinite(W) ||
            error("designer produced no finite certified point")
        return (alpha_best = alpha, pep_best = W,
                n_solves = trace.nsolves, status = "OK")
    catch err
        return (alpha_best = NaN, pep_best = NaN, n_solves = 0,
                status = "ERROR: " * first(sprint(showerror, err), 120))
    end
end

# ── crash-safe Arrow storage ─────────────────────────────────────────────────

"""
    ResultWriter(RowType, dir; batch = 64)

Thread-safe, append-only result store. Rows are buffered and flushed as
self-contained Arrow shards `part-XXXXX.arrow` inside `dir`, so a run that is
killed keeps every already-flushed batch. Resumes without collisions when the
directory already contains shards.
"""
mutable struct ResultWriter{R}
    const dir::String
    const batch::Int
    const lock::ReentrantLock
    buf::Vector{R}
    part::Int
    nrows::Int
end

function ResultWriter(::Type{R}, dir::AbstractString; batch::Int = 64) where {R}
    mkpath(dir)
    part = 0
    for f in readdir(dir)
        mm = match(r"^part-(\d+)\.arrow$", f)
        mm === nothing || (part = max(part, parse(Int, mm.captures[1])))
    end
    return ResultWriter{R}(String(dir), batch, ReentrantLock(), R[], part, 0)
end

"""Write the buffered rows to a new Arrow shard; returns the number written."""
function flush_writer!(w::ResultWriter{R}) where {R}
    local rows::Vector{R}, idx::Int
    lock(w.lock)
    try
        isempty(w.buf) && return 0
        rows = copy(w.buf)
        empty!(w.buf)
        w.part += 1
        idx = w.part
    finally
        unlock(w.lock)
    end
    tmp = joinpath(w.dir, @sprintf("part-%05d.arrow.tmp", idx))
    Arrow.write(tmp, rows)
    mv(tmp, joinpath(w.dir, @sprintf("part-%05d.arrow", idx)); force = true)
    return length(rows)
end

"""Buffer one row, flushing automatically once the batch size is reached."""
function record!(w::ResultWriter{R}, row::R) where {R}
    n = lock(w.lock) do
        push!(w.buf, row)
        w.nrows += 1
        length(w.buf)
    end
    n >= w.batch && flush_writer!(w)
    return nothing
end

close_writer!(w::ResultWriter) = (flush_writer!(w); w.nrows)

"""Concatenate every Arrow shard of `dir` into a single `DataFrame`."""
function load_results(dir::AbstractString)
    isdir(dir) || return DataFrame()
    files = sort(joinpath.(dir,
        filter(f -> occursin(r"^part-\d+\.arrow$", f), readdir(dir))))
    frames = DataFrame[]
    for f in files
        try
            # Own the columns instead of retaining Arrow's memory map. This
            # lets resumed runs and Windows tests close/delete shard files
            # deterministically after loading them.
            push!(frames, DataFrame(Arrow.Table(f); copycols = true))
        catch err
            @warn "skipping unreadable Arrow shard" file = f err
        end
    end
    isempty(frames) && return DataFrame()
    return reduce(vcat, frames)
end

"""
    load_array_results(dir)

Load Arrow shards from a distributed run. Array workers write below
`dir/tasks/task-XXXXX-of-XXXXX`; a local run writes directly in `dir`.
Only shard files are read, so a stale `merged.arrow` is never folded back in.
"""
function load_array_results(dir::AbstractString)
    taskroot = joinpath(dir, "tasks")
    isdir(taskroot) || return load_results(dir)
    frames = DataFrame[]
    for taskdir in sort(filter(isdir, readdir(taskroot; join = true)))
        df = load_results(taskdir)
        nrow(df) == 0 || push!(frames, df)
    end
    isempty(frames) && return DataFrame()
    return reduce(vcat, frames)
end

"""Merge all shards of `dir` into `dir/merged.arrow` (offline-friendly)."""
function merge_parts(dir::AbstractString; out::AbstractString = joinpath(dir, "merged.arrow"))
    df = load_results(dir)
    nrow(df) == 0 && return nothing
    Arrow.write(out, df)
    return out
end

"""Merge every task-local Arrow shard of a Slurm-array run."""
function merge_array_parts(dir::AbstractString;
                           out::AbstractString = joinpath(dir, "merged.arrow"))
    df = load_array_results(dir)
    nrow(df) == 0 && return nothing
    tmp = out * ".tmp"
    Arrow.write(tmp, df)
    mv(tmp, out; force = true)
    return out
end

_key_value(v) = v isa AbstractString ? String(v) : v

"""Set of `cols`-tuples already present in `dir` (used to resume a run)."""
function completed_keys(dir::AbstractString, cols)
    df = load_results(dir)
    nrow(df) == 0 && return Set{Tuple}()
    missingcols = [c for c in cols if !hasproperty(df, c)]
    isempty(missingcols) || error("stored results lack key columns $missingcols")
    return Set{Tuple}(Tuple(_key_value(r[c]) for c in cols) for r in eachrow(df))
end

"""Dump the experiment configuration next to its results."""
function write_config(dir::AbstractString, cfg)
    mkpath(dir)
    open(joinpath(dir, "config.txt"), "w") do io
        println(io, "# written ", Dates.now())
        println(io, "# julia threads = ", Threads.nthreads())
        println(io, repr(merge(cfg, (solver_backend = backend_name(default_backend()),))))
    end
    return nothing
end

# ── parallel helpers ─────────────────────────────────────────────────────────

function _arg_value(args, flag::AbstractString, default)
    i = findfirst(==(flag), args)
    i === nothing && return default
    i < length(args) || error("$flag requires a value")
    return args[i + 1]
end

"""
    array_context([args]) -> (task_id, task_count)

Return a zero-based Slurm-array identity. Command-line `--task-id` and
`--task-count` override `SLURM_ARRAY_TASK_ID` and `PEP_ARRAY_TASK_COUNT`.
Outside Slurm the default is the single task `(0, 1)`.
"""
function array_context(args = ARGS)
    id = parse(Int, _arg_value(args, "--task-id",
                               get(ENV, "SLURM_ARRAY_TASK_ID", "0")))
    count = parse(Int, _arg_value(args, "--task-count",
                                  get(ENV, "PEP_ARRAY_TASK_COUNT", "1")))
    count > 0 || error("task count must be positive, got $count")
    0 <= id < count || error("task id must be in 0:$(count - 1), got $id")
    return (task_id = id, task_count = count)
end

"""Deterministically select the items owned by a zero-based array task."""
function array_partition(items, task_id::Int, task_count::Int)
    0 <= task_id < task_count || error("invalid array task $task_id/$task_count")
    return [item for (i, item) in enumerate(items) if mod(i - 1, task_count) == task_id]
end

"""Collision-free output directory for one array worker."""
function array_task_dir(root::AbstractString, task_id::Int, task_count::Int)
    task_count == 1 && return String(root)
    name = @sprintf("task-%05d-of-%05d", task_id, task_count)
    return joinpath(root, "tasks", name)
end

const COMMON_CLUSTER_CFG_FIELDS = (
    :T_lo, :T_hi,
    :A_lo, :A_hi, :n_A,
    :alpha_lo, :alpha_hi, :n_alpha,
    :kappa_lo, :kappa_hi, :n_kappa,
    :norm_B, :R_x, :ratio_radius, :L_f,
    :facial_reduction, :fr_rtol, :fr_atol, :cert_tol, :gap_tol, :fd_rel_step,
    :batch_rows,
)

_parse_cluster_value(default::Symbol, raw::AbstractString) =
    Symbol(lowercase(strip(raw)))
_parse_cluster_value(default, raw::AbstractString) = parse(typeof(default), raw)

"""
    cluster_overrides(cfg)

Override the common scientific fields of `cfg` from `PEP_*` environment
variables when `PEP_CLUSTER_CONFIG=1`. With that switch absent, return `cfg`
unchanged, which preserves each driver's local defaults.
"""
function cluster_overrides(cfg::NamedTuple)
    get(ENV, "PEP_CLUSTER_CONFIG", "0") == "1" || return cfg
    values = map(keys(cfg)) do key
        default = getproperty(cfg, key)
        key in COMMON_CLUSTER_CFG_FIELDS || return default
        envkey = "PEP_" * uppercase(String(key))
        raw = get(ENV, envkey, nothing)
        raw === nothing ? default : _parse_cluster_value(default, raw)
    end
    return NamedTuple{keys(cfg)}(values)
end

"""
    cluster_optimizers(default)

Read `PEP_OPTIMIZERS` when cluster configuration is enabled. Accepted values
are `som`, `fom`, or a comma-separated combination such as `som,fom`.
"""
function cluster_optimizers(default)
    get(ENV, "PEP_CLUSTER_CONFIG", "0") == "1" || return default
    raw = get(ENV, "PEP_OPTIMIZERS", nothing)
    raw === nothing && return default
    optimizers = unique(Symbol(lowercase(strip(value)))
                        for value in split(raw, ',') if !isempty(strip(value)))
    isempty(optimizers) && error("PEP_OPTIMIZERS must contain som and/or fom")
    invalid = filter(opt -> opt ∉ (:som, :fom), optimizers)
    isempty(invalid) || error("invalid PEP_OPTIMIZERS value(s): $(join(invalid, ", "))")
    return collect(optimizers)
end

"""
    cluster_horizons(cfg)

Return the exact horizon grid. Cluster submissions may provide a
comma-separated `PEP_T_VALUES` list; local runs and older cluster commands
without that variable retain the inclusive `T_lo:T_hi` interval.
"""
function cluster_horizons(cfg::NamedTuple)
    raw = get(ENV, "PEP_CLUSTER_CONFIG", "0") == "1" ?
          strip(get(ENV, "PEP_T_VALUES", "")) : ""
    isempty(raw) && return collect(cfg.T_lo:cfg.T_hi)
    values = Int[]
    for token in split(raw, ',')
        value = strip(token)
        isempty(value) && error("PEP_T_VALUES contains an empty entry")
        push!(values, parse(Int, value))
    end
    all(>(0), values) || error("PEP_T_VALUES must contain positive integers")
    length(unique(values)) == length(values) ||
        error("PEP_T_VALUES must not contain duplicates")
    sort!(values)
    return values
end

"""Validate bounds and sample counts shared by both control drivers."""
function validate_common_cfg(cfg::NamedTuple)
    cfg.T_lo >= 0 || error("T_lo must be nonnegative")
    cfg.T_lo <= cfg.T_hi || error("T_lo must not exceed T_hi")
    cfg.A_lo >= 0 || error("A_lo must be nonnegative")
    cfg.A_lo <= cfg.A_hi || error("A_lo must not exceed A_hi")
    cfg.alpha_lo <= cfg.alpha_hi || error("alpha_lo must not exceed alpha_hi")
    cfg.kappa_lo >= 0 || error("kappa_lo must be nonnegative")
    cfg.kappa_lo <= cfg.kappa_hi || error("kappa_lo must not exceed kappa_hi")
    cfg.kappa_hi <= 1 ||
        error("kappa_hi must not exceed 1 because kappa = mu_f / L_f")
    cfg.n_A > 0 || error("n_A must be positive")
    cfg.n_alpha > 0 || error("n_alpha must be positive")
    cfg.n_kappa > 0 || error("n_kappa must be positive")
    cfg.norm_B >= 0 || error("norm_B must be nonnegative")
    cfg.R_x >= 0 || error("R_x must be nonnegative")
    cfg.ratio_radius >= 0 || error("ratio_radius must be nonnegative")
    cfg.L_f > 0 || error("L_f must be positive")
    cfg.batch_rows > 0 || error("batch_rows must be positive")
    if hasproperty(cfg, :facial_reduction)
        cfg.facial_reduction in (:none, :explicit) ||
            error("facial_reduction must be :none or :explicit")
        cfg.fr_rtol >= 0 || error("fr_rtol must be nonnegative")
        cfg.fr_atol >= 0 || error("fr_atol must be nonnegative")
        cfg.cert_tol >= 0 || error("cert_tol must be nonnegative")
        cfg.gap_tol >= 0 || error("gap_tol must be nonnegative")
        if hasproperty(cfg, :fd_rel_step)
            cfg.fd_rel_step >= 0 || error("fd_rel_step must be nonnegative")
        end
    end
    return cfg
end

"""
    best_sweep_alphas(results) -> Dict

Return the alpha with the smallest finite PEP value for every
`(T, norm_A, kappa)` group in a sweep result table. `results` may be a
DataFrame, an Arrow file, or a run directory. If the directory has no
`merged.arrow`, its available task shards are read (even if the array is
incomplete). Missing groups are intentionally absent so callers can apply a
fallback.
"""
function best_sweep_alphas(df::AbstractDataFrame)
    required = (:T, :norm_A, :kappa, :alpha, :pep)
    missingcols = filter(col -> !hasproperty(df, col), required)
    isempty(missingcols) || error("sweep results lack columns $missingcols")
    best = Dict{Tuple{Int, Float64, Float64}, Tuple{Float64, Float64}}()
    accepted_status = ("OPTIMAL", "SLOW_PROGRESS", "ALMOST_OPTIMAL", "OK")
    for row in eachrow(df)
        hasproperty(df, :status) && String(row.status) ∉ accepted_status && continue
        alpha = Float64(row.alpha)
        pep = Float64(row.pep)
        isfinite(alpha) && isfinite(pep) || continue
        key = (Int(row.T), Float64(row.norm_A), Float64(row.kappa))
        previous = get(best, key, (NaN, Inf))
        pep < previous[2] && (best[key] = (alpha, pep))
    end
    result = Dict{Tuple{Int, Float64, Float64}, Float64}()
    for (key, value) in best
        result[key] = value[1]
    end
    return result
end

function best_sweep_alphas(path::AbstractString)
    if isdir(path)
        file = joinpath(path, "merged.arrow")
        if isfile(file)
            return best_sweep_alphas(DataFrame(Arrow.Table(file); copycols = true))
        end
        shards = load_array_results(path)
        nrow(shards) == 0 &&
            return Dict{Tuple{Int, Float64, Float64}, Float64}()
        return best_sweep_alphas(shards)
    end
    isfile(path) || return Dict{Tuple{Int, Float64, Float64}, Float64}()
    return best_sweep_alphas(DataFrame(Arrow.Table(path); copycols = true))
end

"""Short stable hash of a configuration NamedTuple (run directory tag)."""
cfg_hash(cfg) = string(hash(repr(cfg)) % UInt32; base = 16, pad = 8)

"""Logarithmically spaced grid of `n` points between `lo` and `hi`."""
log_grid(lo::Real, hi::Real, n::Int) =
    n == 1 ? [Float64(lo)] : exp10.(range(log10(lo), log10(hi); length = n))

"""Linearly spaced normalized strong-convexity ratios, including zero."""
kappa_grid(lo::Real, hi::Real, n::Int) =
    n == 1 ? [Float64(lo)] : collect(range(Float64(lo), Float64(hi); length = n))

"""
    pforeach(f, items; nthreads = Threads.nthreads())

Dynamically scheduled `f(item)` over `items`: workers pull the next index from
a shared atomic counter, which keeps threads busy despite the very uneven cost
of the instances (the SDP size grows with the horizon).
"""
function pforeach(f, items; nthreads::Int = Threads.nthreads())
    n = length(items)
    n == 0 && return nothing
    counter = Threads.Atomic{Int}(1)
    @sync for _ in 1:max(1, min(nthreads, n))
        Threads.@spawn begin
            while true
                i = Threads.atomic_add!(counter, 1)
                i > n && break
                f(items[i])
            end
        end
    end
    return nothing
end

end # module
