# # include functions from helper_functions.jl
# # https://docs.julialang.org/en/v1/manual/code-loading/
# include("helper_functions.jl")

# # initialization dict for openloop at nominal conditions
# # params, taken from Zhang et al 2001
# γ = 200.0                     # 1) 1/h, maximum fission rate, when cell m ≥ m_d
# A = sqrt(25/π)*1.0e11         # 2) 1/g, partition function exponential prefactors (i.e., p(m,m'))
# S_l = 0.1                     # 3) g/L, low substrate limit
# K_t = 0.01e-11                # 4) g/g/L, proportionality constant for change in transition mass with substrate levels
# m_t0 = 6.0e-11                # 5) g, initial transition mass (at low substrate levels)
# m_max = 12.0e-11              # 6) g, maximal cell mass allowed
# Y = 0.4                       # 7) g/g, yield coefficient for cell mass to mols substrate
# K_m = 25.0                    # 8) g/L, Monod constant for single cell growth rate
# D = 0.40                      # 9) 1/h, dilution rate in continuous bioreactor operation
# ϵ = 5.0e22                    # 10) g^{-2}, inverse variance for fission rate when cells are in fissioning state
# β = 100.0e22                  # 11) g^{-2}, another inverse variance for partition function
# S_h = 2.0                     # 12) g/L, high substrate limit
# K_d = 2.0e-11                 # 13) g/g/L, proportionality constant for change in critical fissioning mass with substrate levels
# m_d0 = 11.0e-11               # 14) g, initial critical fissioning mass (at low substrate levels)
# m_a = 1.0e-11                 # 15) g, additional mass that mother cell must gain to start fissioning
# μ_m = 5.0e-10                 # 16) g/h, single cell growth rate in limit of high substrate levels
# α = 20.0                      # 17) 1/h, rate constant for delayed adjustment of cell metabolism to substrate levels
# S_f = 25.0                    # 18) g/L, feed concentration of substrate to continuous bioreactor
# N_m = 150                     # 19) -, number of grid points in mass domain

# # IC conditions
# μ_0 = 3.0e-11                 # g, initial cell mass distribution mean
# σ_0 = 1.0e-11                 # g, initial cell mass distribution standard deviation
# N_00 = 1.0e4                  # -, initial total number of cells
# m_edges = gen_m_edges(N_m,0.0,m_max) # collect(range(0.0,m_max,N_m+1))
# mc =  gen_m_centers(m_edges) # ((m_edges[2] - m_edges[1]) / 2.0) .+ m_edges[1:end-1]
# N_IC = N_00 .* N_pdf(mc,μ_0,σ_0)
# S_0 = S_f
# S_prime0 = S_0
# m0 = N_00
# u_0 = [
#     S_0,                      # initial substrate concentration [g/L]
#     S_prime0,                 # initial filtered substrate concentration [g/L]
#     m0                        # from N_IC exact integration [1/L]
# ]

# m_t_test = py"m_td"(S_prime0, m_t0, K_t, S_l, S_h)
# m_d_test = py"m_td"(S_prime0, m_d0, K_d, S_l, S_h)
# Γ0 =py"Γ"(mc,γ,ϵ,m_t_test,m_d_test,m_a)
# Γ0m0 = py"midpt_trapz"((mc[2] - mc[1]) .* ones(length(mc)), (N_IC .* Γ0))

# p_PBE = [
#     γ,
#     A,
#     S_l,
#     K_t,
#     m_t0,
#     m_max,
#     Y,
#     K_m,
#     D,
#     ϵ,
#     β,
#     S_h,
#     K_d,
#     m_d0,
#     m_a,
#     μ_m,
#     α,
#     S_f,
#     N_m,
#     N_IC,
#     Γ0m0
# ]

# initialization_dict_nominal = Dict{String,Any}(
#     "u_0"=>u_0,
#     "N_mt"=>N_IC,
#     "D"=>p_PBE[9],
#     "Sf"=>p_PBE[18],
#     "params"=>p_PBE
# )

# # write to revised jld file
# jld_fp_figureplots = pwd() * "/initialization_dict_nominal.jld"

# # write
# jldopen(jld_fp_figureplots, "w") do file
#     write(
#         file, 
#         "initialization_dict_nominal",
#         initialization_dict_nominal
#         )
# end

# println("Script end.")

println("Script start.")

# open contents of Project.toml (explicitly added packages) and Manifest.toml (recursive dependency packages)
import Pkg; Pkg.activate(".")

# include functions from helper_functions.jl
# https://docs.julialang.org/en/v1/manual/code-loading/
include("helper_functions.jl")


# span the desired range of inputs (D,Sf) creating oscillations, possibly
init_openloop_list = []
D_list = 0.01 .* collect(25:35)
Sf_list = 1.0 .* collect(10:40)

for D in D_list
    for S_f in Sf_list

        # initialization dict for openloop at nominal conditions
        # params, taken from Zhang et al 2001
        γ = 200.0                     # 1) 1/h, maximum fission rate, when cell m ≥ m_d
        A = sqrt(25/π)*1.0e11         # 2) 1/g, partition function exponential prefactors (i.e., p(m,m'))
        S_l = 0.1                     # 3) g/L, low substrate limit
        K_t = 0.01e-11                # 4) g/g/L, proportionality constant for change in transition mass with substrate levels
        m_t0 = 6.0e-11                # 5) g, initial transition mass (at low substrate levels)
        m_max = 12.0e-11              # 6) g, maximal cell mass allowed
        Y = 0.4                       # 7) g/g, yield coefficient for cell mass to mols substrate
        K_m = 25.0                    # 8) g/L, Monod constant for single cell growth rate
        # D = 0.40                      # 9) 1/h, dilution rate in continuous bioreactor operation
        ϵ = 5.0e22                    # 10) g^{-2}, inverse variance for fission rate when cells are in fissioning state
        β = 100.0e22                  # 11) g^{-2}, another inverse variance for partition function
        S_h = 2.0                     # 12) g/L, high substrate limit
        K_d = 2.0e-11                 # 13) g/g/L, proportionality constant for change in critical fissioning mass with substrate levels
        m_d0 = 11.0e-11               # 14) g, initial critical fissioning mass (at low substrate levels)
        m_a = 1.0e-11                 # 15) g, additional mass that mother cell must gain to start fissioning
        μ_m = 5.0e-10                 # 16) g/h, single cell growth rate in limit of high substrate levels
        α = 20.0                      # 17) 1/h, rate constant for delayed adjustment of cell metabolism to substrate levels
        # S_f = 25.0                    # 18) g/L, feed concentration of substrate to continuous bioreactor
        N_m = 150                     # 19) -, number of grid points in mass domain

        # IC conditions
        μ_0 = 3.0e-11                 # g, initial cell mass distribution mean
        σ_0 = 1.0e-11                 # g, initial cell mass distribution standard deviation
        N_00 = 1.0e4                  # -, initial total number of cells
        m_edges = gen_m_edges(N_m,0.0,m_max) # collect(range(0.0,m_max,N_m+1))
        mc =  gen_m_centers(m_edges) # ((m_edges[2] - m_edges[1]) / 2.0) .+ m_edges[1:end-1]
        N_IC = N_00 .* N_pdf(mc,μ_0,σ_0)
        S_0 = S_f
        S_prime0 = S_0
        m0 = N_00
        u_0 = [
            S_0,                      # initial substrate concentration [g/L]
            S_prime0,                 # initial filtered substrate concentration [g/L]
            m0                        # from N_IC exact integration [1/L]
        ]

        m_t_test = py"m_td"(S_prime0, m_t0, K_t, S_l, S_h)
        m_d_test = py"m_td"(S_prime0, m_d0, K_d, S_l, S_h)
        Γ0 =py"Γ"(mc,γ,ϵ,m_t_test,m_d_test,m_a)
        Γ0m0 = py"midpt_trapz"((mc[2] - mc[1]) .* ones(length(mc)), (N_IC .* Γ0))

        p_PBE = [
            γ,
            A,
            S_l,
            K_t,
            m_t0,
            m_max,
            Y,
            K_m,
            D,
            ϵ,
            β,
            S_h,
            K_d,
            m_d0,
            m_a,
            μ_m,
            α,
            S_f,
            N_m,
            N_IC,
            Γ0m0
        ]


        initialization_dict_nominal = Dict{String,Any}(
            "u_0"=>u_0,
            "N_mt"=>N_IC,
            "D"=>p_PBE[9],
            "Sf"=>p_PBE[18],
            "params"=>p_PBE
        )

        push!(init_openloop_list, initialization_dict_nominal)
    end
end

# write to revised jld file
jld_fp_figureplots = pwd() * "/init_openloop_list.jld"

# write
jldopen(jld_fp_figureplots, "w") do file
    write(
        file, 
        "init_openloop_list",
        init_openloop_list
        )
end

println("Script end.")