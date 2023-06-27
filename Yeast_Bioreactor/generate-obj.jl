# include functions from helper_functions.jl
# https://docs.julialang.org/en/v1/manual/code-loading/
include("helper_functions.jl")

# load closed loop results to obtain objective evals
jld_fp_figureplots = pwd() * "/closedloop_DS_SISO_Kc_varied.jld"
closedloop_DS_SISO_Kc_varied = jldopen(jld_fp_figureplots, "r") do file
    read(file, "closedloop_DS_SISO_Kc_varied")
end
# display(closedloop_DS_SISO_Kc_varied)

jld_fp_figureplots = pwd() * "/openloop_nominal_response.jld"
openloop_nominal_response = jldopen(jld_fp_figureplots, "r") do file
    read(file, "openloop_nominal_response")
end

# obtain oscillation half-width
init_SSvals_dict_closedloop = gen_SS_vals(
    openloop_nominal_response["state_values"][1]
)
# println(init_SSvals_dict_closedloop)

# compute objectives but no plotting
obj_res = compute_G_v2(
    closedloop_DS_SISO_Kc_varied;
    previous_p = plot(),
    m0_norm = init_SSvals_dict_closedloop["m0_SS"],
    S_norm = init_SSvals_dict_closedloop["S_SS"],
    m0_weight = 1.0,
    Sf_weight = 1.0,
    γ_1 = 1.0,
    γ_2 = 1.0,
    SISO_label = "(D,S)",
    shade_col = palette(:tab10)[1],
    G_type = "economic",
    plot_type = "none",
    plot_xvar = "Kc",
    filter_yvar_ind = 1,
    reshape_xdim = length(closedloop_DS_SISO_Kc_varied["controller_params"]),
    reshape_ydim = 1
)
# display(obj_res)

println("Script end.")