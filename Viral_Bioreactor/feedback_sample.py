from scipy.integrate import solve_ivp, simpson
import numpy as np
import matplotlib.pyplot as plt

#Plot formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize


"""
Define model
"""

def controlled_model(t,x,u,params,control_plan,control_params):
    #Unpack x, u, params
    #Dictionaries for indexing
    state_assignment = {
        "T": 0,
        "Id": 1,
        "Is": 2,
        "Ic": 3,
        "Vs": 4,
        "Vd": 5
    }

    T = x[0]
    Id = x[1]
    Is = x[2]
    Ic = x[3]
    Vs = x[4]
    Vd  = x[5]

    kc, r = control_params

    
    if control_plan[1] == "Tin":
        Tin =  kc*x[state_assignment[control_plan[0]]] + r
        if Tin < 0.0:
            Tin = 0.0
        elif Tin > 1e7:
            Tin = 1e7
        else:
            Tin = Tin
        D = u[1]
    elif control_plan[1] == "D":
        D =  kc*x[state_assignment[control_plan[0]]] + r
        if D < 0.0:
            D = 0.0
        elif D > 0.1:
            Tin = 1e7
        else:
            D = D
        Tin = u[0]
    else:
        raise RuntimeError("Wrong CV is Specified in control_plan")

    mu = params[0]
    k1 = params[1]
    k2 = params[2]
    k3 = params[3]
    k33 = params[4]
    k4 = params[5]
    f = params[6]

    #Define model
    dTdt = mu*T - k1*(Vs+Vd)*T + D*(Tin-T)
    dIddt = k1*Vd*T - (k1*Vs - mu)*Id - (D)*Id 
    dIsdt = k1*Vs*T -(k1*Vd + k2)*Is -(D)*Is 
    dIcdt = k1*(Vs*Id + Vd*Is) -k2*Ic - (D)*Ic 
    dVsdt = k3*Is - (k1*(T+Id+Is+Ic) + k4 + (D))*Vs 
    dVddt = k33*Ic + f*k3*Is - (k1*(T+Id+Is+Ic) + k4 + (D))*Vd 

    print(t)

    return [dTdt,dIddt,dIsdt,dIcdt,dVsdt,dVddt]

def post_calculator(sol,u,control_plan,control_params):
    #Extract t and state vals from solution
    t_vals = sol.t
    x = [sol.y[0], sol.y[1],sol.y[2],sol.y[3],sol.y[4],sol.y[5]]
    kc,r = control_params
    #Evaluate u
    state_assignment = {
        "T": 0,
        "Id": 1,
        "Is": 2,
        "Ic": 3,
        "Vs": 4,
        "Vd": 5
    }

    if control_plan[1] == "Tin":
        Tin = np.zeros(len(t_vals))
        for i in range(len(t_vals)):
            Tin_ = kc*x[state_assignment[control_plan[0]]][i] + r
            if Tin_ < 0:
                Tin[i] = 0.0
            elif Tin_ > 1e7:
                Tin[i] = 1e7
            else:
                Tin[i] = Tin_
        D = u[1]*np.ones(len(t_vals))
    elif control_plan[1] == "D":
        D = np.zeros(len(t_vals))
        for i in range(len(t_vals)):
            D_ = kc*x[state_assignment[control_plan[0]]][i] + r
            if D_ < 0:
                D[i] = 0.0
            elif D_ > 0.1:
                D[i] = 0.1
            else:
                D[i] = D_
        Tin = u[0]*np.ones(len(t_vals))
    else:
        raise RuntimeError("Wrong CV is Specified in control_plan")
    
    return t_vals, x, Tin, D


"""
Tin and T
"""

model_params = [0.027,2.12e-9,7.13e-3,168,168,0.035,1e-3]
nominal_u = [3e6,0.0396]
x0 = [5e6,0.0,0.0,0.0,1.25e5,0.0]
control_plan = ["T", "Tin"]
t_vec = [0,600]

# params1 = [-100,5e6]
# params2 = [-100, 1.5e7]
# params3 = [-20,5e6]
# params4 = [-20, 1.5e7]




# #Plot case 1
# sol1 = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan,params1),t_vec,x0)
# t_vals1, x1, Tin1, D1 = post_calculator(sol1,nominal_u,control_plan,params1)
# fig1,ax1 = plt.subplots(1,3, figsize=(8,3), num=1)
# ax1[0].semilogy(t_vals1,x1[0],label=r"$T$")
# ax1[0].semilogy(t_vals1,x1[1],label=r"$I_{d}$")
# ax1[0].semilogy(t_vals1,x1[2],label=r"$I_{s}$")
# ax1[0].semilogy(t_vals1,x1[3],label=r"$I_{c}$")
# ax1[0].semilogy(t_vals1,x1[4],label=r"$V_{s}$")
# ax1[0].semilogy(t_vals1,x1[5],label=r"$V_{d}$")
# ax1[1].plot(t_vals1,Tin1)
# ax1[2].plot(t_vals1,D1)
# ax1[0].set_ylabel(r"Concentration")
# ax1[1].set_ylabel(r"$T_{in}$ (cells/mL)")
# ax1[2].set_ylabel(r"$D$ (1/hr)")
# ax1[0].set_xlabel(r"Time (hr)")
# ax1[1].set_xlabel(r"Time (hr)")
# ax1[2].set_xlabel(r"Time (hr)")
# ax1[0].set_ylim(1e-5,1e11)
# ax1[0].legend(loc="lower right", labelspacing=0.2,ncol=3, columnspacing=0.5, handlelength=0.5)
# fig1.tight_layout()
# plt.savefig("kc_100_u0_5e6.png",dpi=300)

# #Plot case 2
# sol2 = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan,params2),t_vec,x0)
# t_vals2, x2, Tin2, D2 = post_calculator(sol2,nominal_u,control_plan,params2)
# fig2,ax2 = plt.subplots(1,3, figsize=(8,3), num=2)
# ax2[0].semilogy(t_vals2,x2[0],label=r"$T$")
# ax2[0].semilogy(t_vals2,x2[1],label=r"$I_{d}$")
# ax2[0].semilogy(t_vals2,x2[2],label=r"$I_{s}$")
# ax2[0].semilogy(t_vals2,x2[3],label=r"$I_{c}$")
# ax2[0].semilogy(t_vals2,x2[4],label=r"$V_{s}$")
# ax2[0].semilogy(t_vals2,x2[5],label=r"$V_{d}$")
# ax2[1].plot(t_vals2,Tin2)
# ax2[2].plot(t_vals2,D2)
# ax2[0].set_ylabel(r"Concentration")
# ax2[1].set_ylabel(r"$T_{in}$ (cells/mL)")
# ax2[2].set_ylabel(r"$D$ (1/hr)")
# ax2[0].set_xlabel(r"Time (hr)")
# ax2[1].set_xlabel(r"Time (hr)")
# ax2[2].set_xlabel(r"Time (hr)")
# ax2[0].set_ylim(1e-5,1e11)
# ax2[0].legend(loc="lower right", labelspacing=0.2,ncol=3, columnspacing=0.5, handlelength=0.5)
# fig2.tight_layout()
# plt.savefig("kc_100_u0_15e7.png",dpi=300)


# #Plot case 3
# sol3 = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan,params3),t_vec,x0)
# t_vals3, x3, Tin3, D3 = post_calculator(sol3,nominal_u,control_plan,params3)
# fig3,ax3 = plt.subplots(1,3, figsize=(8,3), num=3)
# ax3[0].semilogy(t_vals3,x3[0],label=r"$T$")
# ax3[0].semilogy(t_vals3,x3[1],label=r"$I_{d}$")
# ax3[0].semilogy(t_vals3,x3[2],label=r"$I_{s}$")
# ax3[0].semilogy(t_vals3,x3[3],label=r"$I_{c}$")
# ax3[0].semilogy(t_vals3,x3[4],label=r"$V_{s}$")
# ax3[0].semilogy(t_vals3,x3[5],label=r"$V_{d}$")
# ax3[1].plot(t_vals3,Tin3)
# ax3[2].plot(t_vals3,D3)
# ax3[0].set_ylabel(r"Concentration")
# ax3[1].set_ylabel(r"$T_{in}$ (cells/mL)")
# ax3[2].set_ylabel(r"$D$ (1/hr)")
# ax3[0].set_xlabel(r"Time (hr)")
# ax3[1].set_xlabel(r"Time (hr)")
# ax3[2].set_xlabel(r"Time (hr)")
# ax3[0].set_ylim(1e-5,1e11)
# ax3[0].legend(loc="lower right", labelspacing=0.2,ncol=3, columnspacing=0.5, handlelength=0.5)
# fig3.tight_layout()
# plt.savefig("kc_20_u0_5e6.png",dpi=300)

# #Plot case 3
# sol4 = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan,params4),t_vec,x0)
# t_vals4, x4, Tin4, D4 = post_calculator(sol4,nominal_u,control_plan,params4)
# fig4,ax4 = plt.subplots(1,3, figsize=(8,3), num=4)
# ax4[0].semilogy(t_vals4,x4[0],label=r"$T$")
# ax4[0].semilogy(t_vals4,x4[1],label=r"$I_{d}$")
# ax4[0].semilogy(t_vals4,x4[2],label=r"$I_{s}$")
# ax4[0].semilogy(t_vals4,x4[3],label=r"$I_{c}$")
# ax4[0].semilogy(t_vals4,x4[4],label=r"$V_{s}$")
# ax4[0].semilogy(t_vals4,x4[5],label=r"$V_{d}$")
# ax4[1].plot(t_vals4,Tin4)
# ax4[2].plot(t_vals4,D4)
# ax4[0].set_ylabel(r"Concentration")
# ax4[1].set_ylabel(r"$T_{in}$ (cells/mL)")
# ax4[2].set_ylabel(r"$D$ (1/hr)")
# ax4[0].set_xlabel(r"Time (hr)")
# ax4[1].set_xlabel(r"Time (hr)")
# ax4[2].set_xlabel(r"Time (hr)")
# ax4[0].set_ylim(1e-5,1e11)
# ax4[0].legend(loc="lower right", labelspacing=0.2,ncol=3, columnspacing=0.5, handlelength=0.5)
# fig4.tight_layout()
# plt.savefig("kc_20_u0_15e7.png",dpi=300)



"""
T and D
"""

control_plan1 = ["T", "D"]

params1_ = [-2e-8,0.1]
params2_ = [-2e-8, 0.05]
params3_ = [-1e-6,0.25]
params4_ = [-1e-6,0.05]




#Plot case 1
sol1 = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan1,params1_),t_vec,x0)
t_vals1, x1, Tin1, D1 = post_calculator(sol1,nominal_u,control_plan1,params1_)
fig1,ax1 = plt.subplots(1,3, figsize=(8,3), num=1)
ax1[0].semilogy(t_vals1,x1[0],label=r"$T$")
ax1[0].semilogy(t_vals1,x1[1],label=r"$I_{d}$")
ax1[0].semilogy(t_vals1,x1[2],label=r"$I_{s}$")
ax1[0].semilogy(t_vals1,x1[3],label=r"$I_{c}$")
ax1[0].semilogy(t_vals1,x1[4],label=r"$V_{s}$")
ax1[0].semilogy(t_vals1,x1[5],label=r"$V_{d}$")
ax1[1].plot(t_vals1,Tin1)
ax1[2].plot(t_vals1,D1)
ax1[0].set_ylabel(r"Concentration")
ax1[1].set_ylabel(r"$T_{in}$ (cells/mL)")
ax1[2].set_ylabel(r"$D$ (1/hr)")
ax1[0].set_xlabel(r"Time (hr)")
ax1[1].set_xlabel(r"Time (hr)")
ax1[2].set_xlabel(r"Time (hr)")
ax1[0].set_ylim(1e-5,1e11)
ax1[0].legend(loc="lower right", labelspacing=0.2,ncol=3, columnspacing=0.5, handlelength=0.5)
fig1.tight_layout()
plt.savefig("TD_kc_-2e-8_u0_01.png",dpi=300)

#Plot case 2
sol2 = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan1,params2_),t_vec,x0)
t_vals2, x2, Tin2, D2 = post_calculator(sol2,nominal_u,control_plan1,params2_)
fig2,ax2 = plt.subplots(1,3, figsize=(8,3), num=2)
ax2[0].semilogy(t_vals2,x2[0],label=r"$T$")
ax2[0].semilogy(t_vals2,x2[1],label=r"$I_{d}$")
ax2[0].semilogy(t_vals2,x2[2],label=r"$I_{s}$")
ax2[0].semilogy(t_vals2,x2[3],label=r"$I_{c}$")
ax2[0].semilogy(t_vals2,x2[4],label=r"$V_{s}$")
ax2[0].semilogy(t_vals2,x2[5],label=r"$V_{d}$")
ax2[1].plot(t_vals2,Tin2)
ax2[2].plot(t_vals2,D2)
ax2[0].set_ylabel(r"Concentration")
ax2[1].set_ylabel(r"$T_{in}$ (cells/mL)")
ax2[2].set_ylabel(r"$D$ (1/hr)")
ax2[0].set_xlabel(r"Time (hr)")
ax2[1].set_xlabel(r"Time (hr)")
ax2[2].set_xlabel(r"Time (hr)")
ax2[0].set_ylim(1e-5,1e11)
ax2[0].legend(loc="lower right", labelspacing=0.2,ncol=3, columnspacing=0.5, handlelength=0.5)
fig2.tight_layout()
plt.savefig("TD_kc_-2e-8_u0_005.png",dpi=300)



#Plot case 3
sol3 = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan1,params3_),t_vec,x0)
t_vals3, x3, Tin3, D3 = post_calculator(sol3,nominal_u,control_plan1,params3_)
fig3,ax3 = plt.subplots(1,3, figsize=(8,3), num=3)
ax3[0].semilogy(t_vals3,x3[0],label=r"$T$")
ax3[0].semilogy(t_vals3,x3[1],label=r"$I_{d}$")
ax3[0].semilogy(t_vals3,x3[2],label=r"$I_{s}$")
ax3[0].semilogy(t_vals3,x3[3],label=r"$I_{c}$")
ax3[0].semilogy(t_vals3,x3[4],label=r"$V_{s}$")
ax3[0].semilogy(t_vals3,x3[5],label=r"$V_{d}$")
ax3[1].plot(t_vals3,Tin3)
ax3[2].plot(t_vals3,D3)
ax3[0].set_ylabel(r"Concentration")
ax3[1].set_ylabel(r"$T_{in}$ (cells/mL)")
ax3[2].set_ylabel(r"$D$ (1/hr)")
ax3[0].set_xlabel(r"Time (hr)")
ax3[1].set_xlabel(r"Time (hr)")
ax3[2].set_xlabel(r"Time (hr)")
ax3[0].set_ylim(1e-5,1e11)
ax3[0].legend(loc="lower right", labelspacing=0.2,ncol=3, columnspacing=0.5, handlelength=0.5)
fig3.tight_layout()
plt.savefig("TD_kc_-1e-6_u0_25.png",dpi=300)


#Plot case 3
sol4 = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan1,params4_),t_vec,x0)
t_vals4, x4, Tin4, D4 = post_calculator(sol4,nominal_u,control_plan1,params4_)
fig4,ax4 = plt.subplots(1,3, figsize=(8,3), num=4)
ax4[0].semilogy(t_vals4,x4[0],label=r"$T$")
ax4[0].semilogy(t_vals4,x4[1],label=r"$I_{d}$")
ax4[0].semilogy(t_vals4,x4[2],label=r"$I_{s}$")
ax4[0].semilogy(t_vals4,x4[3],label=r"$I_{c}$")
ax4[0].semilogy(t_vals4,x4[4],label=r"$V_{s}$")
ax4[0].semilogy(t_vals4,x4[5],label=r"$V_{d}$")
ax4[1].plot(t_vals4,Tin4)
ax4[2].plot(t_vals4,D4)
ax4[0].set_ylabel(r"Concentration")
ax4[1].set_ylabel(r"$T_{in}$ (cells/mL)")
ax4[2].set_ylabel(r"$D$ (1/hr)")
ax4[0].set_xlabel(r"Time (hr)")
ax4[1].set_xlabel(r"Time (hr)")
ax4[2].set_xlabel(r"Time (hr)")
ax4[0].set_ylim(1e-5,1e12)
ax4[0].legend(loc="lower right", labelspacing=0.2,ncol=3, columnspacing=0.5, handlelength=0.5)
fig4.tight_layout()
plt.savefig("TD_kc_-1e-6_u0_005.png",dpi=300)


plt.show()
