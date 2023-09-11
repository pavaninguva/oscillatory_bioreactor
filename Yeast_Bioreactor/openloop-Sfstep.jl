# include functions from helper_functions.jl
# https://docs.julialang.org/en/v1/manual/code-loading/
include("helper_functions.jl")

# load nominal results
jld_fp_figureplots = pwd() * "/openloop_nominal_response.jld"
openloop_nominal_response = jldopen(jld_fp_figureplots, "r") do file
    read(file, "openloop_nominal_response")
end

# specific_ind depends on the number of entries in openloop_nominal_response 
# (if nominal D=0.4 1/hr and simulated full parameter sweep D∈[0.0,0.5] 1/hr, then at ind=41)
SPECIFIC_IND = 41
initialization_dict_closedloop = extract_init_cond(
    openloop_nominal_response,
    SPECIFIC_IND
)

# do the same with step changes in Sf as in D
TF_CLOSEDLOOP = 20.0
Sf0 = initialization_dict_closedloop["Sf"]
Sf_lo = 0.0
Sf_hi = 200.0
println("Simulating step responses and filters for substrate feed Sf=$(Sf0) g/L.")

# initialize res_dict for open_loop step responses
openloop_Sfnominal_stepresponses = Dict{String, Vector}(
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
        results_dict = openloop_Sfnominal_stepresponses,
        SISO_loop_choice = "none",
        tf = TF_CLOSEDLOOP,
        K_c = 0.0,
        c0 = Sf0,
        u_lo = Sf_lo,
        u_hi = Sf_hi,
        t_control_switch = TF_CLOSEDLOOP,
        ϵ_rel = 0.0
    )

# for stepping, take K_c = 0.0 and c0_new = c0 + stepSf, here stepSf = 20, c0 = 25
stepSf = 15.0
res_stepup = closed_loop_yeastPBM_v2(
        initialization_dict_closedloop;
        results_dict = openloop_Sfnominal_stepresponses,
        SISO_loop_choice = "none",
        tf = TF_CLOSEDLOOP,
        K_c = 0.0,
        c0 = Sf0+stepSf,
        u_lo = Sf_lo,
        u_hi = Sf_hi,
        t_control_switch = TF_CLOSEDLOOP,
        ϵ_rel = 0.0
    )

# for stepping, take K_c = 0.0 and c0_new = c0 + stepD, here stepD = 0.1, c0 = 0.4
res_stepdown = closed_loop_yeastPBM_v2(
        initialization_dict_closedloop;
        results_dict = openloop_Sfnominal_stepresponses,
        SISO_loop_choice = "none",
        tf = TF_CLOSEDLOOP,
        K_c = 0.0,
        c0 = Sf0-stepSf,
        u_lo = Sf_lo,
        u_hi = Sf_hi,
        t_control_switch = TF_CLOSEDLOOP,
        ϵ_rel = 0.0
    )

println("Done")

# write to revised jld file
jld_fp_figureplots = pwd() * "/openloop_Sf_step_response.jld"
openloop_Sf_step_response = openloop_Sfnominal_stepresponses

# write
jldopen(jld_fp_figureplots, "w") do file
    write(
        file, 
        "openloop_Sf_step_response",
        openloop_Sf_step_response
        )
end

# # read
# openloop_Sf_step_response = jldopen(jld_fp_figureplots, "r") do file
#     read(file, "openloop_Sf_step_response")
# end
# display(openloop_Sf_step_response)

println("Script end.")