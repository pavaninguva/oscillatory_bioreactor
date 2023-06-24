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
Define model and post processor
"""

def controlled_model(t,x,u,params,control_plan,kc,r):
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

    return [dTdt,dIddt,dIsdt,dIcdt,dVsdt,dVddt]

def obj_calculator(sol,u,control_plan, kc,r,plot=False):
    #Extract t and state vals from solution
    t_vals = sol.t
    x = [sol.y[0], sol.y[1],sol.y[2],sol.y[3],sol.y[4],sol.y[5]]

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
    
    obj = simpson(np.array(D)*(x[4]+x[5]),t_vals)

    if plot == False:
        return obj
    else:
        return t_vals,x,Tin,D


"""
Run Tin and T
"""


# # Run simulation
model_params = [0.027,2.12e-9,7.13e-3,168,168,0.035,1e-3]
nominal_u = [3e6,0.0396]
x0 = [5e6,0.0,0.0,0.0,1.25e5,0.0]
control_plan = ["T", "Tin"]
kc = -100
r = 5e6
t_vec = [0,400]


# sol1 = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan,kc,r),t_vec,x0)
# t_vals1,x1,Tin1,D1 = obj_calculator(sol1,nominal_u,control_plan,kc, r, plot=True)

# fig1, ax1 = plt.subplots(1,3,figsize=(8,3),num=1)
# ax1[0].semilogy(t_vals1,x1[0],label=r"$T$")
# ax1[0].semilogy(t_vals1,x1[1],label=r"$I_{d}$")
# ax1[0].semilogy(t_vals1,x1[2],label=r"$I_{s}$")
# ax1[0].semilogy(t_vals1,x1[3],label=r"$I_{c}$")
# ax1[0].semilogy(t_vals1,x1[4],label=r"$V_{s}$")
# ax1[0].semilogy(t_vals1,x1[5],label=r"$V_{d}$")
# #Inputs
# ax1[1].plot(t_vals1,Tin1)
# ax1[2].plot(t_vals1,D1)
# #Formatting
# ax1[0].set_ylabel(r"Concentration")
# ax1[1].set_ylabel(r"$T_{in}$ (cells/mL)")
# ax1[2].set_ylabel(r"$D$ (1/hr)")
# ax1[0].set_xlabel(r"Time (hr)")
# ax1[1].set_xlabel(r"Time (hr)")
# ax1[2].set_xlabel(r"Time (hr)")
# ax1[0].set_ylim([1e-5,5e10])
# ax1[0].legend(loc="lower right", labelspacing=0.2,ncol=3, columnspacing=0.5, handlelength=0.5)
# fig1.tight_layout()
# plt.savefig("TinT_feedback_eg.png",dpi=300)


# # #Run contour for Tin and T
# # r_vals = np.linspace(1e6,2e7, 20)
# # kc_vals = np.linspace(-100,-5,20)

# # xx,yy = np.meshgrid(kc_vals,r_vals)
# # Tin_T_opt_vals = np.zeros_like(xx)

# # counter = 0
# # for idx, f in np.ndenumerate(Tin_T_opt_vals):
# #     kc_ = xx[idx]
# #     r_ = yy[idx]
# #     sol_ = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan,kc_,r_),t_vec,x0)
# #     Tin_T_opt_vals[idx] = obj_calculator(sol_,nominal_u,control_plan,kc_, r_, plot=False)
# #     counter = counter + 1
# #     print(counter)

# # fig2 = plt.figure(num=2,figsize=(5,4))
# # plt.pcolormesh(xx,yy,Tin_T_opt_vals,cmap="jet",shading="gouraud")
# # plt.colorbar(label="Objective Function")
# # plt.xlabel(r"$K_{c}$")
# # plt.ylabel(r"$u_{0}$")
# # # plt.yscale("symlog")
# # plt.tight_layout()
# # plt.savefig("TinT_feedback_contour.png",dpi=300)


# """
# Run D and T
# """
control_plan1 = ["T", "D"]
# kc1 = -1e-7
# r1 = 0.1

# sol2 = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan1,kc1,r1),t_vec,x0)
# t_vals2,x2,Tin2,D2 = obj_calculator(sol2,nominal_u,control_plan1,kc1, r1, plot=True)

# fig3, ax3 = plt.subplots(1,3,figsize=(8,3),num=3)
# ax3[0].semilogy(t_vals2,x2[0],label=r"$T$")
# ax3[0].semilogy(t_vals2,x2[1],label=r"$I_{d}$")
# ax3[0].semilogy(t_vals2,x2[2],label=r"$I_{s}$")
# ax3[0].semilogy(t_vals2,x2[3],label=r"$I_{c}$")
# ax3[0].semilogy(t_vals2,x2[4],label=r"$V_{s}$")
# ax3[0].semilogy(t_vals2,x2[5],label=r"$V_{d}$")
# #Inputs
# ax3[1].plot(t_vals2,Tin2)
# ax3[2].plot(t_vals2,D2)
# #Formatting
# ax3[0].set_ylabel(r"Concentration")
# ax3[1].set_ylabel(r"$T_{in}$ (cells/mL)")
# ax3[2].set_ylabel(r"$D$ (1/hr)")
# ax3[0].set_xlabel(r"Time (hr)")
# ax3[1].set_xlabel(r"Time (hr)")
# ax3[2].set_xlabel(r"Time (hr)")
# ax3[0].set_ylim([1e-5,5e10])
# ax3[0].legend(loc="lower right", labelspacing=0.2,ncol=3, columnspacing=0.5, handlelength=0.5)
# fig3.tight_layout()
# plt.savefig("DT_feedback_eg.png",dpi=300)


kc_vals = -np.logspace(-6,-8,25)
r_vals = np.linspace(0.25,0.03,25)

xx,yy = np.meshgrid(kc_vals,r_vals)
Tin_D_opt_vals = np.zeros_like(xx)

counter = 0
for idx, f in np.ndenumerate(Tin_D_opt_vals):
    kc_ = xx[idx]
    r_ = yy[idx]
    sol_ = solve_ivp(lambda t,x: controlled_model(t,x,nominal_u,model_params,control_plan1,kc_,r_),t_vec,x0)
    Tin_D_opt_vals[idx] = obj_calculator(sol_,nominal_u,control_plan1,kc_, r_, plot=False)
    counter = counter + 1
    print(counter)

fig4 = plt.figure(num=4,figsize=(5,4))
plt.pcolormesh(xx,yy,Tin_D_opt_vals,cmap="jet",shading="gouraud")
plt.colorbar(label="Objective Function")
plt.xlabel(r"$K_{c}$")
plt.xscale("symlog")
plt.ylabel(r"$u_{0}$")
plt.tight_layout()
plt.savefig("DT_feedback_contour.png",dpi=300)



plt.show()
