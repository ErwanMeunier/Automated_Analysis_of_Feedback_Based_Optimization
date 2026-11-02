# Re-run specific archived FOM initializers; use a small case by default.
import Pkg
Pkg.activate(joinpath(@__DIR__,"..");io=devnull)
include(joinpath(@__DIR__,"..","src","ReviewData.jl"))
include(joinpath(@__DIR__,"..","src","PEPControl.jl"))
using .ReviewData, .PEPControl, DataFrames, Arrow

function main()
    opts=Dict("T"=>"3","norm-A"=>"0.3","kappa"=>"0","iters"=>"historical",
              "solver"=>get(ENV,"PEP_SOLVER","auto"),
              "out"=>joinpath(ROOT,"results","replay_design.arrow"))
    for i in 1:2:length(ARGS)
        i<length(ARGS) || error("$(ARGS[i]) needs a value")
        key=replace(ARGS[i],r"^--"=>"")
        haskey(opts,key) || error("unknown option $(ARGS[i])")
        opts[key]=ARGS[i+1]
    end
    ENV["PEP_SOLVER"]=opts["solver"]
    design=load_reference(corrections=nothing).design
    horizons=parse.(Int,split(opts["T"],","))
    mask=in.(design.T,Ref(horizons))
    opts["norm-A"]=="all" || (mask .&= isapprox.(design.norm_A,parse(Float64,opts["norm-A"]);atol=1e-10))
    opts["kappa"]=="both" || (mask .&= design.kappa .== parse(Float64,opts["kappa"]))
    selected=design[mask,:]
    isempty(selected) && error("no archived groups match the selection")
    result=selected[1:0,:]
    special=Set(((0.0,0.0),(0.1,0.0),(0.2,0.0),(0.0,0.1),(0.1,0.1),(0.4,0.1)))
    for original in eachrow(selected)
        niter=opts["iters"]=="historical" ?
            (original.T==50 && (round(original.norm_A;digits=1),original.kappa) in special ? 60 : 30) :
            parse(Int,opts["iters"])
        inst=ControlInstance(T=Int(original.T),norm_A=Float64(original.norm_A),
            norm_B=Float64(original.norm_B),R_x=Float64(original.R_x),
            R_u=Float64(original.R_u),mu_f=Float64(original.mu_f),L_f=Float64(original.L_f))
        started=time()
        r=best_positive_constant_alpha(compile_control(inst),nsteps(inst),original.alpha_init;
            iters=niter,facial_reduction=:explicit,cert_tol=3e-6,gap_tol=1e-5)
        row=merge(NamedTuple(original),(alpha_best=r.alpha_best,pep_best=r.pep_best,
            n_solves=r.n_solves,status=r.status,solve_time=time()-started))
        push!(result,row;cols=:setequal)
        @info "archived design replay" T=inst.T norm_A=inst.norm_A kappa=original.kappa iterations=niter status=r.status
    end
    out=abspath(opts["out"]);mkpath(dirname(out))
    Arrow.write(out,result)
    println("Wrote ",nrow(result)," replay row(s) to ",relpath(out,ROOT))
end
main()
