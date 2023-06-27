# include functions from helper_functions.jl
# https://docs.julialang.org/en/v1/manual/code-loading/
include("helper_functions.jl")

# load initialization_dict
jld_fp_figureplots = pwd() * "/initialization_dict_nominal.jld"
initialization_dict_nominal = jldopen(jld_fp_figureplots, "r") do file
    read(file, "initialization_dict_nominal")
end

# simulate at nominal conditions
TF_NOMINAL = 100.0
D0 = initialization_dict_nominal["D"]
D_lo = 0.0
D_hi = 1.0
println("Simulating nominal conditions for dilution rate D=$(D0) 1/h.")

# initialize res_dict for open_loop step responses
openloop_nominal = Dict{String, Vector}(
            "cell_mass_distributions"=>[],
            "state_values"=>[],
            "solution_wall_times"=>[],
            "D_values"=>[],
            "Sf_values"=>[],
            "params"=>[],
            "solution_timepts"=>[],
            "controller_params"=>[]
        )

res_nominal = closed_loop_yeastPBM_v2(
        initialization_dict_nominal;
        results_dict = openloop_nominal,
        SISO_loop_choice = "none",
        tf = TF_NOMINAL,
        K_c = 0.0,
        c0 = D0,
        u_lo = D_lo,
        u_hi = D_hi,
        ϵ_rel = 0.0
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

# # read
# openloop_D_step_response = jldopen(jld_fp_figureplots, "r") do file
#     read(file, "openloop_D_step_response")
# end
# display(openloop_D_step_response)

println("Script end.")