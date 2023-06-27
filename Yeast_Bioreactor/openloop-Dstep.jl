# include functions from helper_functions.jl
# https://docs.julialang.org/en/v1/manual/code-loading/
include("helper_functions.jl")

# load nominal results
jld_fp_figureplots = pwd() * "/openloop_nominal_response.jld"
openloop_nominal_response = jldopen(jld_fp_figureplots, "r") do file
    read(file, "openloop_nominal_response")
end

# specific_ind depends on the number of entries in
SPECIFIC_IND = 1
initialization_dict_closedloop = extract_init_cond(
    openloop_nominal_response,
    SPECIFIC_IND
)

# check if wrapped function v2 is working for step changes in D
TF_CLOSEDLOOP = 20.0
D0 = initialization_dict_closedloop["D"]
D_lo = 0.0
D_hi = 1.0
println("Simulating step responses and filters for dilution rate D=$(D0) 1/h.")

# initialize res_dict for open_loop step responses
openloop_Dnominal_stepresponses = Dict{String, Vector}(
            "cell_mass_distributions"=>[],
            "state_values"=>[],
            "solution_wall_times"=>[],
            "D_values"=>[],
            "Sf_values"=>[],
            "params"=>[],
            "solution_timepts"=>[],
            "controller_params"=>[]
        )

res_no_step = closed_loop_yeastPBM_v2(
        initialization_dict_closedloop;
        results_dict = openloop_Dnominal_stepresponses,
        SISO_loop_choice = "none",
        tf = TF_CLOSEDLOOP,
        K_c = 0.0,
        c0 = D0,
        u_lo = D_lo,
        u_hi = D_hi,
        t_control_switch = TF_CLOSEDLOOP,
        ϵ_rel = 0.0
    )

# for stepping, take K_c = 0.0 and c0_new = c0 + stepD, here stepD = 0.1, c0 = 0.4
stepD = 0.1
res_stepup = closed_loop_yeastPBM_v2(
        initialization_dict_closedloop;
        results_dict = openloop_Dnominal_stepresponses,
        SISO_loop_choice = "none",
        tf = TF_CLOSEDLOOP,
        K_c = 0.0,
        c0 = D0+stepD,
        u_lo = D_lo,
        u_hi = D_hi,
        t_control_switch = TF_CLOSEDLOOP,
        ϵ_rel = 0.0
    )

# for stepping, take K_c = 0.0 and c0_new = c0 + stepD, here stepD = 0.1, c0 = 0.4
res_stepdown = closed_loop_yeastPBM_v2(
        initialization_dict_closedloop;
        results_dict = openloop_Dnominal_stepresponses,
        SISO_loop_choice = "none",
        tf = TF_CLOSEDLOOP,
        K_c = 0.0,
        c0 = D0-stepD,
        u_lo = D_lo,
        u_hi = D_hi,
        t_control_switch = TF_CLOSEDLOOP,
        ϵ_rel = 0.0
    )

println("Done")

# write to revised jld file
jld_fp_figureplots = pwd() * "/openloop_D_step_response.jld"
openloop_D_step_response = openloop_Dnominal_stepresponses

# write
jldopen(jld_fp_figureplots, "w") do file
    write(
        file, 
        "openloop_D_step_response",
        openloop_D_step_response
        )
end

# # read
# openloop_D_step_response = jldopen(jld_fp_figureplots, "r") do file
#     read(file, "openloop_D_step_response")
# end
# display(openloop_D_step_response)

println("Script end.")