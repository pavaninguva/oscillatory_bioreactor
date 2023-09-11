# open contents of Project.toml (explicitly added packages) and Manifest.toml (recursive dependency packages)
import Pkg; Pkg.activate(".")

using Distributed
addprocs(72, topology=:master_worker, exeflags="--project=$(Base.active_project())"); # modify to desired number of cores for parallel computing

# include functions from helper_functions.jl
# https://docs.julialang.org/en/v1/manual/code-loading/
@everywhere include("helper_functions.jl")

# load list of initialization_dicts
jld_fp_figureplots = pwd() * "/init_openloop_list.jld"
init_openloop_list = jldopen(jld_fp_figureplots, "r") do file
    read(file, "init_openloop_list")
end

# simulate at nominal conditions using parallel map
TF_NOMINAL = 200.0 
D_lo = 0.0
D_hi = 0.5
println("Simulating in parallel nominal conditions for dilution rate and substrate feed concentrations.")

res_nominal = pmap(
    (init_dict) -> closed_loop_yeastPBM_v2(
        init_dict;
        SISO_loop_choice = "none",
        tf = TF_NOMINAL,
        K_c = 0.0,
        c0 = init_dict["D"],
        u_lo = D_lo,
        u_hi = D_hi,
        ϵ_rel = 0.0
    ),
    init_openloop_list
)

println("Done")

# write to revised jld file
jld_fp_figureplots = pwd() * "/openloop_nominal_response.jld"
openloop_nominal_response = res_nominal

# write
jldopen(jld_fp_figureplots, "w") do file
    write(
        file, 
        "openloop_nominal_response",
        openloop_nominal_response
        )
end

println("Script end.")