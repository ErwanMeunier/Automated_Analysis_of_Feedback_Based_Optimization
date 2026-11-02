import Pkg
const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT; io=devnull)
using TOML, Dates

function options(args)
    opts = Dict("config"=>joinpath(ROOT,"config","smoke.toml"),
                "out"=>joinpath(ROOT,"results","smoke"), "stage"=>"all",
                "solver"=>get(ENV,"PEP_SOLVER","auto"))
    for i in 1:2:length(args)
        i < length(args) || error("$(args[i]) needs a value")
        key = replace(args[i], r"^--"=>"")
        haskey(opts,key) || error("unknown option $(args[i])")
        opts[key] = args[i+1]
    end
    return opts
end

function main()
    opts = options(ARGS)
    cfg = TOML.parsefile(opts["config"])
    horizons = Int.(cfg["horizons"])
    !isempty(horizons) && all(>(0), horizons) && allunique(horizons) ||
        error("horizons must be distinct positive integers")
    kappas = Float64.(cfg["kappas"])
    length(kappas) in (1,2) && all(k -> 0 <= k <= 1, kappas) && allunique(kappas) ||
        error("provide one kappa, or two distinct kappa endpoints")
    cfg["A_lo"] >= 0 && cfg["A_hi"] < 1 && cfg["A_lo"] <= cfg["A_hi"] ||
        error("paper reproduction requires 0 <= A_lo <= A_hi < 1")
    cfg["n_A"] > 0 && cfg["n_alpha"] > 0 || error("sample counts must be positive")
    opts["solver"] in ("auto","mosek","clarabel") || error("invalid solver")
    stage = opts["stage"]
    ngroups = length(horizons) * Int(cfg["n_A"]) * length(kappas)
    stage == "count" && (println(ngroups); return)
    stage in ("all","validate","sweep","design","analytic","merge-sweep","merge-design","merge-analytic") ||
        error("unknown stage $stage")
    out = abspath(opts["out"])
    mkpath(out)
    # One shared snapshot also guards separate-stage invocations from mixing
    # a new design configuration with an old sweep in the same directory.
    plan=Dict("scientific_configuration"=>cfg,"requested_solver"=>opts["solver"])
    plan_path=joinpath(out,"experiment.toml")
    if !isfile(plan_path)
        temporary,io=mktemp(out)
        try
            TOML.print(io,plan); close(io)
            try
                mv(temporary,plan_path;force=false)
            catch
                isfile(plan_path) || rethrow()
            end
        finally
            isopen(io) && close(io)
            isfile(temporary) && rm(temporary)
        end
    end
    TOML.parsefile(plan_path)==plan ||
        error("output directory has a different configuration/solver; select a new --out")
    stage=="validate" && (println("Scientific configuration validated; groups=$ngroups"); return)
    env = Dict{String,String}("PEP_CLUSTER_CONFIG"=>"1", "PEP_SOLVER"=>opts["solver"],
        "PEP_T_VALUES"=>join(horizons,","), "PEP_OPTIMIZERS"=>"fom",
        "PEP_KAPPA_LO"=>string(minimum(kappas)), "PEP_KAPPA_HI"=>string(maximum(kappas)),
        "PEP_N_KAPPA"=>string(length(kappas)), "PEP_DESIGN_ITERS"=>string(cfg["design_iters"]),
        "PEP_DESIGN_INIT_FROM_SWEEP"=>"1", "PEP_DESIGN_FALLBACK_ALPHA"=>"0.25",
        "PEP_SWEEP_RESULTS"=>joinpath(out,"sweep","merged.arrow"),
        "PEP_BATCH_ROWS"=>"1", "OPENBLAS_NUM_THREADS"=>"1", "GKSwstype"=>"100")
    for key in ("A_lo","A_hi","n_A","alpha_lo","alpha_hi","n_alpha","L_f",
                "norm_B","R_x","ratio_radius","facial_reduction","cert_tol","gap_tol","fd_rel_step")
        env["PEP_"*uppercase(key)] = string(cfg[key])
    end
    # Record only scientific parameters, never the inherited environment or a
    # license/account/host setting. Every array worker has a unique manifest.
    task = get(ENV,"SLURM_ARRAY_TASK_ID",get(ENV,"PEP_ARRAY_TASK_ID","local"))
    manifest = merge(cfg, Dict("requested_solver"=>opts["solver"],
        "julia_version"=>string(VERSION), "stage"=>stage,
        "started_utc"=>string(now(UTC)), "groups"=>ngroups))
    manifest_path=joinpath(out,"run_$(stage)_$(task).toml")
    if isfile(manifest_path)
        previous=TOML.parsefile(manifest_path)
        all(get(previous,key,nothing)==value for (key,value) in cfg) &&
            get(previous,"requested_solver",nothing)==opts["solver"] ||
            error("output directory has a different configuration/solver; select a new --out")
    end
    open(manifest_path,"w") do io
        TOML.print(io,manifest)
    end
    function execute(script, args...)
        cmd = `$(Base.julia_cmd()) --project=$ROOT --threads=$(Threads.nthreads()) $(joinpath(@__DIR__,script)) $args`
        process=run(ignorestatus(addenv(cmd, env...)))
        success(process) || error("$script failed; see the stage's diagnostic output")
    end
    for current in (stage == "all" ? ("sweep","design","analytic") : (stage,))
        if startswith(current,"merge-")
            kind = replace(current,"merge-"=>"")
            execute("merge_array_results.jl",joinpath(out,kind),string(ngroups))
        elseif current == "analytic"
            execute("evaluate_analytic_from_pep.jl",joinpath(out,"design","merged.arrow"),
                "--method","theorem1","--facial-reduction",cfg["facial_reduction"],
                "--cert-tol",string(cfg["cert_tol"]),"--gap-tol",string(cfg["gap_tol"]),
                "--out",joinpath(out,"analytic"))
        else
            script = current == "sweep" ? "sweep_alpha.jl" : "design_best_alpha.jl"
            execute(script,"--out",joinpath(out,current))
        end
    end
end
main()
