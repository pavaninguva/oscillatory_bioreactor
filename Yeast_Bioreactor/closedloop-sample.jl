# include functions from helper_functions.jl
# https://docs.julialang.org/en/v1/manual/code-loading/
include("helper_functions.jl")

# load nominal results
jld_fp_figureplots = pwd() * "/openloop_nominal_response.jld"
openloop_nominal_response = jldopen(jld_fp_figureplots, "r") do file
    read(file, "openloop_nominal_response")
end

# specific_ind depends on the number o# specific_ind depends on the number of entries in openloop_nominal_response 
# (if nominal D=0.4 1/hr and simulated full parameter sweep D∈[0.0,0.5] 1/hr, then at ind=41)f entries in
SPECIFIC_IND = 41
initialization_dict_closedloop = extract_init_cond(
    openloop_nominal_response,
    SPECIFIC_IND
)

# do the same with step changes in Sf as in D
TF_CLOSEDLOOP = 20.0
D0 = initialization_dict_closedloop["D"]
D_lo = 0.0
D_hi = 1.0
println("Simulating closed loop responses for dilution rate D=$(D0) 1/h.")

# initialize res_dict for open_loop step responses
closedloop_DS_SISO_Kc_varied = Dict{String, Vector}(
            "cell_mass_distributions"=>[],
            "state_values"=>[],
            "solution_wall_times"=>[],
            "D_values"=>[],
            "Sf_values"=>[],
            "params"=>[],
            "solution_timepts"=>[],
            "controller_params"=>[]
        )

# span a probable K_c range (K_c = ΔD/ΔS ≈ 0.1 / 0.5 ≈ 0.2)
# washout happens above K_c = 0.2
# gain is positive, with feedback law u=Kc*x+u0, need to take negative of gain to get stable process
for i in collect(1:10)
    res_tuned = closed_loop_yeastPBM_v2(
            initialization_dict_closedloop;
            results_dict = closedloop_DS_SISO_Kc_varied,
            SISO_loop_choice = "(D,S)",
            tf = TF_CLOSEDLOOP,
            K_c = -i*0.2,
            c0 = D0,
            u_lo = D_lo,
            u_hi = D_hi,
            t_control_switch = 0.0,
            ϵ_rel = 0.0
        )
end

println("Done")

# write to revised jld file
jld_fp_figureplots = pwd() * "/closedloop_DS_SISO_Kc_varied.jld"

# write
jldopen(jld_fp_figureplots, "w") do file
    write(
        file, 
        "closedloop_DS_SISO_Kc_varied",
        closedloop_DS_SISO_Kc_varied
        )
end

println("Script end.")