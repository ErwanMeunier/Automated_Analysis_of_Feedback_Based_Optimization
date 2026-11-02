import Pkg
Pkg.activate(joinpath(@__DIR__,"..");io=devnull)
using Test, Arrow, DataFrames, SHA, TOML
include(joinpath(@__DIR__,"..","src","ReviewData.jl"))
using .ReviewData
let i=findfirst(==("--solver"),ARGS)
    if i!==nothing
        i<length(ARGS) || error("--solver needs a value")
        ENV["PEP_SOLVER"]=ARGS[i+1]
    end
end

@testset "Reference data integrity and correction policy" begin
    provenance=TOML.parsefile(joinpath(ROOT,"data","reference","PROVENANCE.toml"))
    for record in provenance["files"]
        path=joinpath(ROOT,"data","reference",record["file"])
        @test bytes2hex(open(sha256,path))==record["sha256"]
        @test nrow(DataFrame(Arrow.Table(path)))==record["rows"]
    end
    tables=load_reference(corrections=nothing)
    @test nrow(tables.sweep)==42000
    @test nrow(tables.design)==440
    @test nrow(tables.analytic)==440
    @test all(tables.design.norm_A .< 1)
    @test Set(tables.design.kappa)==Set((0.0,0.1))
    # The paper's zero-strong-convexity case must never be replaced by the supplement.
    @test all(tables.design.mu_f[tables.design.kappa .== 0].==0)
    original=DataFrame(T=[3],norm_A=[0.3],kappa=[0.0],value=[1.0],status=["OK"])
    failure=DataFrame(T=[3],norm_A=[0.3],kappa=[0.0],value=[99.0],status=["UNCERTIFIED"])
    @test overlay(original,failure).value==[1.0]
    success=DataFrame(T=[3],norm_A=[0.3],kappa=[0.0],value=[0.9],status=["OK"])
    @test overlay(original,success).value==[0.9]
    physical=copy(original); physical.L_f=[1.0]
    wrong_units=copy(success); wrong_units.L_f=[2.0]
    @test_throws ErrorException overlay(physical,wrong_units)
end

if !("--data-only" in ARGS)
    include(joinpath(@__DIR__,"..","src","PEPControl.jl"))
    include(joinpath(@__DIR__,"..","src","PEPWorstCase.jl"))
    using .PEPControl, .PEPWorstCase
    const PD=PEPControl.PEPDesign
    @testset "Small SDP, certificate checks, and physical scaling" begin
        backend=PD.default_backend()
        @info "Testing SDP backend" solver=PD.backend_name(backend)
        for kappa in (0.0,0.1)
            inst=ControlInstance(T=2,norm_A=0.3,mu_f=kappa)
            cp=compile_control(inst)
            sol=PD.solve_pep(cp,fill(0.1,nsteps(inst));backend,
                facial_reduction=:explicit,cert_tol=3e-6,gap_tol=1e-5)
            @test sol.diagnostics.certified
            @test isfinite(sol.obj) && sol.obj>0
            scaled=ControlInstance(T=2,norm_A=0.3,mu_f=2kappa,L_f=2.0)
            scaled_sol=PD.solve_pep(compile_control(scaled),fill(0.05,nsteps(scaled));backend,
                facial_reduction=:explicit,cert_tol=3e-6,gap_tol=1e-5)
            @test scaled_sol.diagnostics.certified
            @test isapprox(scaled_sol.obj,2sol.obj;rtol=2e-4)
            @test isapprox(theorem1_uniform_step_size(0.3,1.0,2.0),
                           theorem1_uniform_step_size(0.3,1.0,1.0)/2;rtol=1e-12)
        end
    end
end
