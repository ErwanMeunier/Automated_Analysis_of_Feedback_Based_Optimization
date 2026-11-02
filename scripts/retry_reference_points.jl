# Bounded repairs. Original reference Arrow files are never modified.
import Pkg
Pkg.activate(joinpath(@__DIR__,".."); io=devnull)
include(joinpath(@__DIR__,"..","src","ReviewData.jl"))
include(joinpath(@__DIR__,"..","src","PEPControl.jl"))
include(joinpath(@__DIR__,"..","src","PEPWorstCase.jl"))
using .ReviewData, .PEPControl, .PEPWorstCase
using Arrow, CSV, DataFrames, Dates
const PD=PEPControl.PEPDesign

function main()
    opts=Dict("max-T"=>"12","limit"=>"2","design-iters"=>"3","seconds"=>"60",
        "solver"=>get(ENV,"PEP_SOLVER","auto"),"out"=>joinpath(ROOT,"data","corrections"))
    for i in 1:2:length(ARGS)
        i<length(ARGS) || error("$(ARGS[i]) needs a value")
        key=replace(ARGS[i],r"^--"=>"")
        haskey(opts,key) || error("unknown option $(ARGS[i])")
        opts[key]=ARGS[i+1]
    end
    maxT,limit,iters=parse.(Int,(opts["max-T"],opts["limit"],opts["design-iters"]))
    min(maxT,limit,iters)>=1 || error("retry limits must be positive")
    ENV["PEP_SOLVER"]=opts["solver"]
    ENV["PEP_SOLVE_TIME_LIMIT"]=opts["seconds"]
    tables=load_reference(corrections=opts["out"])
    audit=audit_reference(tables)
    out=abspath(opts["out"]); mkpath(out)
    history_path=joinpath(out,"retry_attempts.csv")
    history=isfile(history_path) ? CSV.read(history_path,DataFrame;
        types=Dict(:recorded_utc=>String,:message=>String)) : DataFrame()
    attempts=NamedTuple[]
    attempted=0
    for row in eachrow(audit), kind in (:analytic,:design)
        row[Symbol(kind,"_ok")] && continue
        table=getproperty(tables,kind)
        row_index=findfirst(r->point_key(r)==point_key(row),eachrow(table))
        row_index===nothing && error("no parameter row for $(point_key(row))")
        original=table[row_index,:]
        started=time()
        result="DEFERRED"; message="outside bounded retry horizon/count"
        gap=scalar=cone=NaN; value=alpha=NaN; solver="not run"
        if row.T<=maxT && attempted<limit
            attempted+=1
            try
                inst=ControlInstance(T=Int(row.T),norm_A=Float64(row.norm_A),
                    norm_B=Float64(original.norm_B),R_x=Float64(original.R_x),
                    R_u=Float64(original.R_u),mu_f=Float64(original.mu_f),L_f=Float64(original.L_f))
                compiled=compile_control(inst)
                solver=PD.backend_name(PD.default_backend())
                alpha=theorem1_uniform_step_size(inst.norm_A,inst.norm_B,inst.L_f)
                corrected=DataFrame([NamedTuple(original)])
                if kind==:design
                    designed=best_positive_constant_alpha(compiled,nsteps(inst),alpha;
                        iters,facial_reduction=:explicit,cert_tol=3e-6,gap_tol=1e-5)
                    isfinite(designed.alpha_best) || error(designed.status)
                    alpha=designed.alpha_best
                    corrected[1,:alpha_init]=theorem1_uniform_step_size(inst.norm_A,inst.norm_B,inst.L_f)
                    corrected[1,:alpha_best]=alpha
                    corrected[1,:pep_best]=designed.pep_best
                    corrected[1,:n_solves]=designed.n_solves
                end
                solution=PD.solve_pep(compiled,fill(alpha,nsteps(inst));warn=false,
                    facial_reduction=:explicit,cert_tol=3e-6,gap_tol=1e-5)
                diag=solution.diagnostics
                gap,scalar,cone=diag.relative_gap,diag.max_scalar_violation,diag.max_cone_violation
                value=solution.obj
                result=diag.certified && isfinite(value) ? "OK" : "UNCERTIFIED"
                message=string(solution.status)
                if kind==:analytic
                    for (field,v) in ((:alpha_analytic,alpha),(:pep_analytic,value),
                        (:analytic_relative_gap,gap),(:analytic_scalar_violation,scalar),
                        (:analytic_cone_violation,cone),(:analytic_solver_status,message))
                        corrected[1,field]=v
                    end
                else
                    corrected[1,:pep_best]=value
                end
                corrected[1,:status]=result
                corrected[1,:solve_time]=time()-started
                file=kind==:design ? "merged_design.arrow" : "merged_analytic_theorem1.arrow"
                path=joinpath(out,file)
                prior=isfile(path) ? DataFrame(Arrow.Table(path);copycols=true) : corrected[1:0,:]
                prior=prior[[point_key(r)!=point_key(row) for r in eachrow(prior)],:]
                Arrow.write(path,vcat(prior,corrected;cols=:setequal))
            catch err
                result="ERROR"
                message=first(sprint(showerror,err),300)
                message=replace(message,r"[A-Za-z]:[\\/][^\n]*"=>"[local path omitted]")
            end
        end
        push!(attempts,(T=Int(row.T),norm_A=Float64(row.norm_A),kappa=Float64(row.kappa),
            kind=String(kind),solver,status=result,alpha,value,relative_gap=gap,
            scalar_violation=scalar,cone_violation=cone,elapsed_seconds=time()-started,
            message,recorded_utc=string(now(UTC))))
        CSV.write(history_path,vcat(history,DataFrame(attempts);cols=:union))
        @info "reference retry" kind T=row.T norm_A=row.norm_A kappa=row.kappa status=result
    end
    remaining=audit_reference(load_reference(corrections=out))
    CSV.write(joinpath(out,"unresolved.csv"),remaining[remaining.needs_repair,:])
    println("Attempted $attempted point(s); remaining unresolved groups: ",count(remaining.needs_repair))
end
main()
