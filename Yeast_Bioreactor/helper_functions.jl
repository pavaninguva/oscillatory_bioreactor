# open contents of Project.toml (explicitly added packages) and Manifest.toml (recursive dependency packages)
import Pkg; Pkg.activate(".")

# compile the relevant packages
using Printf 
using Conda
using PyCall
using JLD
using LaTeXStrings
using BenchmarkTools
using OrdinaryDiffEq
using DiffEqCallbacks
using Documenter
using Plots; pyplot(html_output_format=:png) # activates `pyplot()` backend instead of `gr()` default
using Random

# wrapper for PyClaw routine
py"""
# import relevant packages
import numpy as np
from scipy.stats import norm
from clawpack import pyclaw
from clawpack import riemann

# define helper functions
N_pdf_py = lambda x,mu,sigma : norm.pdf(x, loc = mu, scale = sigma)
midpt_trapz = lambda delx_list,y_c : np.dot(delx_list,y_c)
k_m = lambda S_prime, μ_m, K_m: (μ_m * S_prime)/(K_m + S_prime)

# define for both transition mass and critical mass
def m_td(S_prime,m_td0,K_td,S_l,S_h):
    if S_prime < S_l:
        return m_td0 + K_td*(S_l - S_h)
    elif (S_prime > S_l) and (S_prime < S_h):
        return m_td0 + K_td*(S_prime - S_h)
    else:
        return m_td0

# define fission rate function
# UPDATE 06/29/23: caught bug in fission rate def'n causing γ spike in Γ results for N_m = 150 at m = (m_t+m_a)
@np.vectorize
def Γ(m,γ,ϵ,m_t,m_d,m_a):
    if m <= m_t+m_a: # used to be m < m_t+m_a, changed to m <= (m_t+m_a)
        return 0.0
    elif (m > m_t+m_a) and (m < m_d):
        return γ*np.exp(-ϵ*(m-m_d)**2)
    else:
        return γ

# partition function
def p_mm(m,m_prime,m_t,m_a,A,β):
    if m > m_prime:
        return 0.0
    elif m_prime < (m_t + m_a):
        return 0.0
    else:
        return A*(np.exp(-β*(m-m_t)**2)+np.exp(-β*(m-m_prime+m_t)**2))

def kernel_mat(m_list, N_m, m_t, m_a, A, β, Γ):

    p_mm_mat = np.zeros((N_m,N_m))

    for j in np.arange(0,N_m):
        for i in np.arange(0,N_m):
            p_mm_mat[j,i] = p_mm(m = m_list[i], m_prime = m_list[j], m_t=m_t, m_a=m_a, A=A, β=β)
    
    # np.multiply is an element-wise multiplication of the two input arrays, i.e., explicit broadcasting
    kernel_m = np.multiply(np.outer(Γ, np.ones(N_m)), p_mm_mat) # mother cell info on columns
    return kernel_m

# for generating the RHS of the solver, but also fed directly to solver
def dq_linhom_src(solver,state,dt):
    
    N_mt = state.q[0,:]
    Γ = state.aux[0,:] # Γ
    mc = state.grid.m.centers
    del_mc_list = np.multiply(np.ones(len(mc)), (mc[1]-mc[0]))

    # compute integral expression
    F = np.zeros(state.problem_data['N_m'])
    for i in range(1,state.problem_data['N_m']+1):
        F[i-1] = midpt_trapz(del_mc_list, np.multiply(state.aux[i,:], N_mt))

    # compute full RHS source term
    dq = np.empty(state.q.shape)
    D = state.problem_data['D']
    dq[0,:] = dt * (-(D + Γ) * N_mt + 2.0 * F)   
    return dq

def step_yeast_PBM_v2(S_prime, p_PBE, N_IC, dt):

    γ = p_PBE[0]                   # 1/hr, maximum fission rate, when cell m ≥ m_d
    A = p_PBE[1]                   # 1/g, partition function exponential prefactors (i.e., p(m,m'))
    S_l = p_PBE[2]                 # g/L, low substrate limit
    K_t = p_PBE[3]                 # g/g/L, proportionality constant for change in transition mass with substrate levels
    m_t0 = p_PBE[4]                # g, initial transition mass (at low substrate levels)
    m_max = p_PBE[5]               # g, maximal cell mass allowed
    Y = p_PBE[6]                   # g/g, yield coefficient for cell mass to mols substrate
    K_m = p_PBE[7]                 # g/L, Monod constant for single cell growth rate
    D = p_PBE[8]                   # 1/h, dilution rate in continuous bioreactor operation
    ϵ = p_PBE[9]                   # g^{-2}, inverse variance for fission rate when cells are in fissioning state
    β = p_PBE[10]                  # g^{-2}, another inverse variance for partition function
    S_h = p_PBE[11]                # g/L, high substrate limit
    K_d = p_PBE[12]                # g/g/L, proportionality constant for change in critical fissioning mass with substrate levels
    m_d0 = p_PBE[13]               # g, initial critical fissioning mass (at low substrate levels)
    m_a = p_PBE[14]                # g, additional mass that mother cell must gain to start fissioning
    μ_m = p_PBE[15]                # g/h, single cell growth rate in limit of high substrate levels
    α = p_PBE[16]                  # 1/h, rate constant for delayed adjustment of cell metabolism to substrate levels
    S_f = p_PBE[17]                # g/L, feed concentration of substrate to continuous bioreactor
    N_m = int(p_PBE[18])           # number of grid points in mass domain
    
    # constant filtered substrate concentration
    m_t_test = m_td(S_prime, m_td0=m_t0, K_td=K_t, S_l=S_l, S_h=S_h)
    m_d_test = m_td(S_prime, m_td0=m_d0, K_td=K_d, S_l=S_l, S_h=S_h)
    m_ta_test = m_t_test + m_a

    # set up solver
    solver = pyclaw.SharpClawSolver1D(riemann.advection_1D_py.advection_1D)
    solver.kernel_language = 'Python'
    solver.weno_order = 5
    solver.time_integrator = 'SSP104'
    solver.dq_src = dq_linhom_src
    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[0]=pyclaw.BC.extrap
    solver.aux_bc_upper[0]=pyclaw.BC.extrap

    # domain
    m = pyclaw.Dimension(0.0,m_max,N_m,name='m')
    domain = pyclaw.Domain(m)
    num_aux=1+N_m
    state = pyclaw.State(domain,solver.num_eqn,num_aux)

    # parameters
    state.problem_data['u'] = k_m(S_prime,μ_m,K_m)  # growth rate
    state.problem_data['D'] = D                     # dilution rate
    state.problem_data['N_m'] = N_m                 # mesh size

    # initialize q and aux vars
    # aux vars are for holding fission rate Γ vector and kernel K(m,m',S') matrix
    mc = state.grid.m.centers
    Γ_res_test = Γ(m = mc, γ=γ, ϵ=ϵ, m_t=m_t_test, m_d=m_d_test, m_a=m_a)
    kernel_test = kernel_mat(m_list=mc, N_m=N_m, m_t=m_t_test, m_a=m_a, A=A, β=β, Γ=Γ_res_test)
    state.q[0,:] = N_IC
    state.q[0,0] = 0.0 # set lowest mass mesh point state value to 0 to help with extrapolating BC on LHS of domain
    state.aux[0,:] = Γ_res_test
    for i in range(1,N_m+1):
        state.aux[i,:] = kernel_test[:,i-1]

    # Controller
    claw = pyclaw.Controller()
    claw.keep_copy = True # saves copy of outputs to claw.frames
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.tfinal = dt
    claw.output_style = 2 # for stepping to exact times desired given out_times
    #claw.num_output_times = 1
    claw.out_times = np.array([dt/2.0,dt]) # output solution only at final step
    claw.verbosity = False
    
    # run and compute zeroth order moment
    claw.run()
    del_mc_list = np.multiply(np.ones(len(claw.centers[0])), (claw.centers[0][1]-claw.centers[0][0]))
    weighted_zero_moment = midpt_trapz(del_mc_list,claw.frames[1].q[0,:] * Γ_res_test)
    return (claw.frames[1].q[0,:], weighted_zero_moment) # grab last timepoint result for N(m,t) and zeroth moment
"""

# ODE approach with computation offloaded to parameter vector
function substrateODEs_v2!(
        du::Vector{Float64},
        u::Vector{Float64},
        p::Vector{Any},
        t::Float64
    )

    # unpack only the necessary params for ODE simulations
    D = p[9]
    S_f = p[18]
    Y = p[7]
    α = p[17]
    μ_m = p[16]
    K_m = p[8]
    Γm0 = p[end]

    # update ODE RHS
    du[1] = D*(S_f - u[1]) - (μ_m*u[2])/(K_m+u[2])/Y * u[3]
    du[2] = α*(u[1] - u[2])
    du[3] = -D * u[3] + Γm0
end

# Gaussian + edge / center creation
N_pdf(x::Union{AbstractFloat,AbstractArray}, mu::AbstractFloat, sigma::AbstractFloat) = @. 1/sigma/sqrt(2*pi) * exp( -(x - mu)^2 / 2 / sigma^2) # normal distribution for h(x) or h(m) depending on whether (x,mu_x,sigma_x) or (m,mu_m,sigma_m) are fed
gen_m_edges(N_m,m_lo,m_hi) = collect(range(m_lo,m_hi,N_m+1))
gen_m_centers(m_edges) = ((m_edges[2] - m_edges[1]) / 2.0) .+ m_edges[1:end-1]

function P_controller(K_c,x,c0,ϵ,u_lo,u_hi)
    u = K_c*(x+ϵ)+c0
    if u >= u_hi
        return u_hi
    elseif u <= u_lo
        return u_lo
    else
        return u
    end
end

function closed_loop_yeastPBM_v2(
        initialization_dict;
        results_dict = Dict{String, Vector}(
            "cell_mass_distributions"=>[],
            "state_values"=>[],
            "solution_wall_times"=>[],
            "D_values"=>[],
            "Sf_values"=>[],
            "params"=>[],
            "solution_timepts"=>[],
            "controller_params"=>[]
        ),
        SISO_loop_choice = "none",
        tf = 10.0,
        K_c = 0.0,
        c0 = 0.0,
        u_lo = 0.0,
        u_hi = Inf,
        t_control_switch = tf,
        ϵ_rel = 0.0
    )
    
    # check if keys of passed dicts are correct
    true_dict_keys = [
        "state_values",
        "cell_mass_distributions",
        "solution_wall_times",
        "D_values",
        "Sf_values",
        "params",
        "solution_timepts",
        "controller_params"
    ]
    
    @assert Set(collect(keys(results_dict))
        ) == Set(true_dict_keys) "Closed loop simulation keys do not match."
    
    true_SISO_options = [
        "none",
        "(D,m0)",
        "(D,S)",
        "(Sf,m0)",
        "(Sf,S)"
    ]
    
    @assert (SISO_loop_choice in true_SISO_options
        ) "Chosen loop $(SISO_loop_choice) not available from $(true_SISO_options)."
    
    # save controller params as dict
    controller_params = Dict{String, Float64}(
        "K_c" => K_c,
        "c0" => c0,
        "t_control_switch" => t_control_switch,
        "ϵ_rel" => ϵ_rel
    )
    
    # unpack params and save params for dynamics and controls
    # deepcopy the parameters vector so that mutation does not propagate to the original initialization dict
    p_PBE = deepcopy(initialization_dict["params"])
    push!(results_dict["params"], p_PBE)
    push!(results_dict["controller_params"], controller_params)

    γ = p_PBE[1]                    # 1) 1/hr, maximum fission rate, when cell m ≥ m_d
    A = p_PBE[2]                    # 2) 1/g, partition function exponential prefactors (i.e., p(m,m'))
    S_l = p_PBE[3]                  # 3) g/L, low substrate limit
    K_t = p_PBE[4]                  # 4) g/g/L, proportionality constant for change in transition mass with substrate levels
    m_t0 = p_PBE[5]                 # 5) g, initial transition mass (at low substrate levels)
    m_max = p_PBE[6]                # 6) g, maximal cell mass allowed
    Y = p_PBE[7]                    # 7) g/g, yield coefficient for cell mass to mols substrate
    K_m = p_PBE[8]                  # 8) g/L, Monod constant for single cell growth rate
    D = p_PBE[9]                    # 9) 1/h, dilution rate in continuous bioreactor operation
    ϵ = p_PBE[10]                   # 10) g^{-2}, inverse variance for fission rate when cells are in fissioning state
    β = p_PBE[11]                   # 11) g^{-2}, another inverse variance for partition function
    S_h = p_PBE[12]                 # 12) g/L, high substrate limit
    K_d = p_PBE[13]                 # 13) g/g/L, proportionality constant for change in critical fissioning mass with substrate levels
    m_d0 = p_PBE[14]                # 14) g, initial critical fissioning mass (at low substrate levels)
    m_a = p_PBE[15]                 # 15) g, additional mass that mother cell must gain to start fissioning
    μ_m = p_PBE[16]                 # 16) g/h, single cell growth rate in limit of high substrate levels
    α = p_PBE[17]                   # 17) 1/h, rate constant for delayed adjustment of cell metabolism to substrate levels
    S_f = p_PBE[18]                 # 18) g/L, feed concentration of substrate to continuous bioreactor
    N_m = p_PBE[19]                 # 19) -, number of grid points in mass domain
    N_IC = p_PBE[20]                # 20) 1/g/L, initial condition for cell number distribution in mass
    Γ0m0 = p_PBE[21]                # 21) 1/h/L, weighted substrate consumption kernel
    
    # IC conditions and numerical params
    u_0 = initialization_dict["u_0"]
    tspan = (0.0,tf)
    
    # pose the problem, no callbacks necessary, no saveat necessary (since controlling Δt)
    prob_openloop = ODEProblem(
        substrateODEs_v2!,
        u_0,
        tspan,
        p_PBE
    )

    integrator = init(
        prob_openloop, 
        Rodas5(autodiff = false),
        reltol = 1e-8,
        abstol = 1e-8
    )
    
    # coupling time scale manually chosen
    # point is to decrease the amount of times at which coupling is forced
    # also cannot push the timestep too low as this will produce memory issues...
    Δt = 1.0e-2
    N_t = length(tspan[1]:Δt:tspan[2])
    times = Vector{Float64}(undef,N_t)
    parsed_N_mt = Matrix{Float64}(undef,N_t,N_m)
    u_list = Vector{Vector{Float64}}(undef,N_t)
    D_closedloop = Vector{Float64}(undef,N_t)
    Sf_closedloop = Vector{Float64}(undef,N_t)
    
    # sequentialism in solve, order of updates is control action, then PBM, then DAE system
    res_time = @elapsed begin 
        for (i,t) in enumerate(tspan[1]:Δt:tspan[2])

            # save some values for plotting, importantly time steps, 
            times[i] = integrator.t
            parsed_N_mt[i,:] = integrator.p[end-1]
            u_list[i] = copy(integrator.u)
            
            # update param vals for correct mesh and control parameters
            ϵ_applied = ϵ_rel * randn() # for adding white noise to m0 sensor
            K_c_switch = K_c * (integrator.t > t_control_switch)
            
            if SISO_loop_choice in ["(D,m0)","(D,S)"]
                
                integrator.p[9] = P_controller(
                    K_c_switch,
                    ifelse(
                        SISO_loop_choice == "(D,m0)",
                        copy(integrator.u[3]),
                        copy(integrator.u[1])
                    ),
                    c0,
                    ϵ_applied,
                    u_lo,
                    u_hi
                )
                
            elseif SISO_loop_choice in ["(Sf,m0)","(Sf,S)"]
                
                integrator.p[18] = P_controller(
                    K_c_switch,
                    ifelse(
                        SISO_loop_choice == "(Sf,m0)",
                        copy(integrator.u[3]),
                        copy(integrator.u[1])
                    ),
                    c0,
                    ϵ_applied,
                    u_lo,
                    u_hi
                )
            
            # don't modify anything if no SISO loop specified
            else
            end
            
            # save both current D or Sf input delivered from controller
            D_closedloop[i] = copy(integrator.p[9])
            Sf_closedloop[i] = copy(integrator.p[18])
            
            N_mt_res, Γm0_res = py"step_yeast_PBM_v2"(
                integrator.u[2],       # S′
                integrator.p[1:end-2], # p_PBE
                integrator.p[end-1],   # N_mt from previous timestamp
                Δt                     # dt for stepping PBM
            )
            
            integrator.p[end-1] = N_mt_res
            integrator.p[end] = Γm0_res
            OrdinaryDiffEq.step!(integrator, Δt, true)

        end
    end
    
    # progress bar
    println("Controller gain=$(K_c) and Controller bias=$(c0) completed in $(res_time)")
    
    # update res_dict
    push!(results_dict["cell_mass_distributions"], parsed_N_mt)
    push!(results_dict["state_values"], mapreduce(permutedims, vcat, u_list)) # u_res
    push!(results_dict["solution_wall_times"], res_time)
    push!(results_dict["D_values"], D_closedloop)
    push!(results_dict["Sf_values"], Sf_closedloop)
    push!(results_dict["solution_timepts"], times)
    
    return results_dict
    
end

function gen_SS_vals(u_res)

    tp_ind = -500
    min_ylim_S = minimum(u_res[end+tp_ind:end,1])
    max_ylim_S = maximum(u_res[end+tp_ind:end,1])
    mid_ylim_S = (max_ylim_S + min_ylim_S) / 2.0

    min_ylim_m0 = minimum(u_res[end+tp_ind:end,3])
    max_ylim_m0 = maximum(u_res[end+tp_ind:end,3])
    mid_ylim_m0 = (max_ylim_m0 + min_ylim_m0) / 2.0

    return Dict(
        "m0_SS"=>mid_ylim_m0,
        "S_SS"=>mid_ylim_S
    )
end

# cumulative trapezoidal rule
function cumtrapz(X::AbstractVector, Y::AbstractVector)
    # Check matching vector length
    @assert length(X) == length(Y)
    
    # Initialize Output
    out = 0.0

    # Iterate over arrays
    @inbounds for i in 2:length(X)
        out += 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
    end
    # Return output
    return out
end

function compute_G_v2(
        results_dict;
        previous_p = plot(),
        m0_norm = 1.0,
        S_norm = 1.0,
        m0_weight = 1.0,
        Sf_weight = 1.0,
        γ_1 = 1.0,
        γ_2 = 1.0,
        SISO_label = "",
        shade_col = "black",
        G_type = "economic",
        plot_type = "none",
        plot_xvar = "Kc",
        filter_yvar_ind = 1,
        reshape_xdim = 1,
        reshape_ydim = 1
    )

    @assert G_type in [
        "economic",
        "controller",
        "setpoint-tracking"
        ] "Requested objective type $(G_type) not available from $(["economic", "controller"])"

    @assert plot_type in ["none", "2d", "3d"] "$(plot_type) neither 2d nor 3d."
    @assert plot_xvar in ["Kc", "c0"] "$(plot_xvar) neither Kc nor c0."

    # check if keys of passed dicts are correct
    true_dict_keys = [
        "state_values",
        "cell_mass_distributions",
        "solution_wall_times",
        "D_values",
        "Sf_values",
        "params",
        "solution_timepts",
        "controller_params"
    ]

    @assert Set(collect(keys(results_dict))
        ) == Set(true_dict_keys) "Closed loop simulation keys do not match."

    # initialize
    G_res1 = []
    G_res2 = []
    G_res3 = []
    K_c_list = []
    c0_list = []

    for (i,res) in enumerate(results_dict["state_values"])
        t_list = results_dict["solution_timepts"][i]
        Sf_norm = results_dict["params"][i][18]
        m0_norm_list = results_dict["state_values"][i][:,3] ./ m0_norm
        Sf_norm_list = results_dict["Sf_values"][i] ./ Sf_norm
        S_norm_list = results_dict["state_values"][i][:,1] ./ S_norm
        
        # first derivative upwind loses one mesh space d.o.f.
        dm0_dt2 = Vector{Float64}(undef,length(t_list))
        dm0_dt2[1] = 0.0
        dm0_dt2[2:end] .= ((m0_norm_list[2:end] .- m0_norm_list[1:end-1])./ (t_list[2:end] .- t_list[1:end-1])).^2

        dS_dt2 = Vector{Float64}(undef,length(t_list))
        dS_dt2[1] = 0.0
        dS_dt2[2:end] .= ((S_norm_list[2:end] .- S_norm_list[1:end-1])./ (t_list[2:end] .- t_list[1:end-1])).^2

        g1 = results_dict["D_values"][i] .* (m0_weight .* m0_norm_list .- Sf_weight .* Sf_norm_list)
        g2 = -γ_1 .* dm0_dt2 .- γ_2 .* dS_dt2
        g3 = -(abs.(m0_norm_list .- 1.0) .+ abs.(S_norm_list .- 1.0))
        
    #         println(cumtrapz(t_list, results_dict["D_values"][i] .* m0_weight .* m0_norm_list), " Biomass contribution")
    #         println(cumtrapz(t_list, results_dict["D_values"][i] .* Sf_weight .* Sf_norm_list), " Sf contribution")        
    #         println(cumtrapz(t_list, γ_1 .* dm0_dt2)," m0 oscillation contribution")
    #         println(cumtrapz(t_list, γ_2 .* dS_dt2)," S oscillation contribution")
    #         println(cumtrapz(t_list, g3)," SP tracking contribution")
        
        push!(G_res1,cumtrapz(t_list,g1))
        push!(G_res2,cumtrapz(t_list,g2))
        push!(G_res3,cumtrapz(t_list,g3))
        push!(K_c_list, results_dict["controller_params"][i]["K_c"])
        push!(c0_list, results_dict["controller_params"][i]["c0"])
        
    end

    # normalize K_c_list by maximal value
    K_c_list_norm = abs.(K_c_list) ./ maximum(abs.(K_c_list))
    c0_list_norm = abs.(c0_list) ./ maximum(abs.(c0_list))

    # check if reshaping is possible, and reshape to do 2d + 3d plotting properly
    @assert (
        reshape_xdim * reshape_ydim) == length(K_c_list_norm
        ) "Requested reshape dims (y_dim,x_dim) = ($(reshape_ydim),$(reshape_xdim)) do not match K_c size $(length(K_c_list_norm))"
    @assert (
        reshape_xdim * reshape_ydim) == length(c0_list_norm
        ) "Requested reshape dims (y_dim,x_dim) = ($(reshape_ydim),$(reshape_xdim)) do not match c0 size $(length(c0_list_norm))"

    # do the reshaping
    K_c_mat_norm = reshape(K_c_list_norm, (reshape_ydim,reshape_xdim))
    c0_mat_norm = reshape(c0_list_norm, (reshape_ydim,reshape_xdim))
    G_res1_mat = reshape(G_res1, (reshape_ydim,reshape_xdim))
    G_res2_mat = reshape(G_res2, (reshape_ydim,reshape_xdim))
    G_res3_mat = reshape(G_res3, (reshape_ydim,reshape_xdim))

    # plot
    if plot_type == "2d"
        if plot_xvar == "Kc"
            if G_type == "economic"

                # draw lines connecting the scatter points
                plot!(
                    previous_p,
                    K_c_mat_norm[filter_yvar_ind,:],
                    G_res1_mat[filter_yvar_ind,:],
                    label = "",
                    color = shade_col,
                    linewidth = 1
                )
                
                scatter!(
                    previous_p,
                    K_c_mat_norm[filter_yvar_ind,:],
                    G_res1_mat[filter_yvar_ind,:],
                    label = SISO_label,
                    xlabel = L"$\tilde{K}_c$",
                    ylabel = L"$G_1(t_f)$",
                    color = shade_col,
                    legend = :bottomleft,
                    framestyle = :box,
                    grid = :off,
                    xlim = (0.0,1.05)
                )
            elseif G_type == "controller"
                
                # draw lines connecting the scatter points
                plot!(
                    previous_p,
                    K_c_mat_norm[filter_yvar_ind,:],
                    G_res2_mat[filter_yvar_ind,:],
                    label = "",
                    color = shade_col,
                    linewidth = 1
                )
                
                scatter!(
                    previous_p,
                    K_c_mat_norm[filter_yvar_ind,:],
                    G_res2_mat[filter_yvar_ind,:],
                    label = SISO_label,
                    xlabel = L"$\tilde{K}_c$",
                    ylabel = L"$G_2(t_f)$",
                    color = shade_col,
                    legend = :bottomleft,
                    framestyle = :box,
                    grid = :off,
                    xlim = (0.0,1.05)
                )
            else
                
                # draw lines connecting the scatter points
                plot!(
                    previous_p,
                    K_c_mat_norm[filter_yvar_ind,:],
                    G_res3_mat[filter_yvar_ind,:],
                    label = "",
                    color = shade_col,
                    linewidth = 1
                )
                
                scatter!(
                    previous_p,
                    K_c_mat_norm[filter_yvar_ind,:],
                    G_res3_mat[filter_yvar_ind,:],
                    label = SISO_label,
                    xlabel = L"$\tilde{K}_c$",
                    ylabel = L"$G_3(t_f)$",
                    color = shade_col,
                    legend = :bottomleft,
                    framestyle = :box,
                    grid = :off,
                    xlim = (0.0,1.05)
                )
            end
            
        elseif plot_xvar == "c0"
            
            if G_type == "economic"

                # draw lines connecting the scatter points
                plot!(
                    previous_p,
                    c0_mat_norm[:,filter_yvar_ind],
                    G_res1_mat[:,filter_yvar_ind],
                    label = "",
                    color = shade_col,
                    linewidth = 1
                )
                
                scatter!(
                    previous_p,
                    c0_mat_norm[:,filter_yvar_ind],
                    G_res1_mat[:,filter_yvar_ind],
                    label = SISO_label,
                    xlabel = L"$\tilde{c}_0$",
                    ylabel = L"$G_1(t_f)$",
                    color = shade_col,
                    legend = :bottomleft,
                    framestyle = :box,
                    grid = :off,
                    xlim = (0.0,1.05)
                )
                
            elseif G_type == "controller"
                
                # draw lines connecting the scatter points
                plot!(
                    previous_p,
                    c0_mat_norm[:,filter_yvar_ind],
                    G_res2_mat[:,filter_yvar_ind],
                    label = "",
                    color = shade_col,
                    linewidth = 1
                )
                
                scatter!(
                    previous_p,
                    c0_mat_norm[:,filter_yvar_ind],
                    G_res2_mat[:,filter_yvar_ind],
                    label = SISO_label,
                    xlabel = L"$\tilde{c}_0$",
                    ylabel = L"$G_2(t_f)$",
                    color = shade_col,
                    legend = :bottomleft,
                    framestyle = :box,
                    grid = :off,
                    xlim = (0.0,1.05)
                )
                
            else
                
                # draw lines connecting the scatter points
                plot!(
                    previous_p,
                    c0_mat_norm[:,filter_yvar_ind],
                    G_res3_mat[:,filter_yvar_ind],
                    label = "",
                    color = shade_col,
                    linewidth = 1
                )
                
                scatter!(
                    previous_p,
                    c0_mat_norm[:,filter_yvar_ind],
                    G_res3_mat[:,filter_yvar_ind],
                    label = SISO_label,
                    xlabel = L"$\tilde{c}_0$",
                    ylabel = L"$G_3(t_f)$",
                    color = shade_col,
                    legend = :bottomleft,
                    framestyle = :box,
                    grid = :off,
                    xlim = (0.0,1.05)
                )
                
            end
        
        end
        return previous_p

    elseif plot_type == "3d"

        if G_type == "economic"
            plot!(
                previous_p,
                K_c_mat_norm,
                c0_mat_norm,
                convert.(Float64, G_res1_mat),
                color = shade_col,
                label = "",
                xlabel = L"$\tilde{K}_c$",
                ylabel = L"$\tilde{c}_0$",
                zlabel = L"$G_1(t_f)$",
                xlim = (0.0,1.05),
                ylim = (0.0,1.05),
                legend = false
            )
        elseif G_type == "controller"
            plot!(
                previous_p,
                K_c_mat_norm,
                c0_mat_norm,
                convert.(Float64, G_res2_mat),
                color = shade_col,
                label = "",
                xlabel = L"$\tilde{K}_c$",
                ylabel = L"$\tilde{c}_0$",
                zlabel = L"$G_2(t_f)$",
                xlim = (0.0,1.05),
                ylim = (0.0,1.05),
                legend = false
            )
        else
            plot!(
                previous_p,
                K_c_mat_norm,
                c0_mat_norm,
                convert.(Float64, G_res3_mat),
                color = shade_col,
                label = "",
                xlabel = L"$\tilde{K}_c$",
                ylabel = L"$\tilde{c}_0$",
                zlabel = L"$G_3(t_f)$",
                xlim = (0.0,1.05),
                ylim = (0.0,1.05),
                legend = false
            )
        end
        return previous_p

    # if "none" for plotting, return the raw results without re-shaping and normalizing
    else
        objective_dict = Dict{String,Vector}(
            "Kc_normalized"=>K_c_list,
            "c0_normalized"=> c0_list,
            "economic"=> G_res1,
            "controller"=> G_res2,
            "setpoint-tracking"=> G_res3
        )
        return objective_dict
    end

end


# results_dict = Dict{String, Vector}(
#             "cell_mass_distributions"=>[],
#             "state_values"=>[],
#             "solution_wall_times"=>[],
#             "D_values"=>[],
#             "Sf_values"=>[],
#             "params"=>[],
#             "solution_timepts"=>[],
#             "controller_params"=>[]
#         ),
# extract initialization conditions for a simulation
# specific_ind refers to a positional entry in the containers of the res_dict added
function extract_init_cond(res_dict_list,specific_ind)

    # create specific endpoint reconstruction dict with appropriate updates to, e.g., Γm0
    p_PBE = res_dict_list[specific_ind]["params"][1]
    γ = p_PBE[1]                   # 1/hr, maximum fission rate, when cell m ≥ m_d
    A = p_PBE[2]                   # 1/g, partition function exponential prefactors (i.e., p(m,m'))
    S_l = p_PBE[3]                 # g/L, low substrate limit
    K_t = p_PBE[4]                 # g/g/L, proportionality constant for change in transition mass with substrate levels
    m_t0 = p_PBE[5]                # g, initial transition mass (at low substrate levels)
    m_max = p_PBE[6]               # g, maximal cell mass allowed
    Y = p_PBE[7]                   # g/g, yield coefficient for cell mass to mols substrate
    K_m = p_PBE[8]                 # g/L, Monod constant for single cell growth rate
    D = p_PBE[9]                   # 1/h, dilution rate in continuous bioreactor operation
    ϵ = p_PBE[10]                  # g^{-2}, inverse variance for fission rate when cells are in fissioning state
    β = p_PBE[11]                  # g^{-2}, another inverse variance for partition function
    S_h = p_PBE[12]                # g/L, high substrate limit
    K_d = p_PBE[13]                # g/g/L, proportionality constant for change in critical fissioning mass with substrate levels
    m_d0 = p_PBE[14]               # g, initial critical fissioning mass (at low substrate levels)
    m_a = p_PBE[15]                # g, additional mass that mother cell must gain to start fissioning
    μ_m = p_PBE[16]                # g/h, single cell growth rate in limit of high substrate levels
    α = p_PBE[17]                  # 1/h, rate constant for delayed adjustment of cell metabolism to substrate levels
    S_f = p_PBE[18]                # g/L, feed concentration of substrate to continuous bioreactor
    N_m = p_PBE[19]                # number of grid points in mass domain
    m_edges = collect(range(0.0,m_max,N_m+1))
    mc = ((m_edges[2] - m_edges[1]) / 2.0) .+ m_edges[1:end-1]
    N_IC = res_dict_list[specific_ind]["cell_mass_distributions"][1][end,:]
    u_res = res_dict_list[specific_ind]["state_values"][1]
    m_t_test = py"m_td"(u_res[end,2], m_t0, K_t, S_l, S_h)
    m_d_test = py"m_td"(u_res[end,2], m_d0, K_d, S_l, S_h)
    Γ0 = py"Γ"(mc,γ,ϵ,m_t_test,m_d_test,m_a)  
    p_PBE[20] = N_IC
    p_PBE[21] = py"midpt_trapz"((mc[2] - mc[1]) .* ones(length(mc)), (p_PBE[20] .* Γ0))

    u_0 = res_dict_list[specific_ind]["state_values"][1][end,:]

    init_dict = Dict{String,Any}(
        "u_0"=>u_0,
        "N_mt"=>N_IC,
        "D"=>p_PBE[9],
        "Sf"=>p_PBE[18],
        "params"=>p_PBE
    )

    return init_dict
end