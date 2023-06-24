from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

#Plot formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size="14")
plt.rc('legend', fontsize=10)
plt.rc("axes", titlesize=12)


def model(t,x,u, control_plan, kc, r):
    X, S, P, Z, W = x

    #State_assignment
    state_assignment = {
        "X":0,
        "S":1,
        "P":2
    }

    if control_plan[1] == "Sf":
        Sf = kc*x[state_assignment[control_plan[0]]] + r
        if Sf < 0.0:
            Sf = 0.0
        elif Sf > 200.0:
            Sf = 200.0
        else:
            Sf = Sf
        D = u[0]
    elif control_plan[1] == "D":
        D = kc*x[state_assignment[control_plan[0]]] + r
        if D < 0.0:
            D = 0.0
        elif D > 0.06:
            D = 0.06
        else:
            D = D
        Sf = u[1]

    print(D)
    print(Sf)

    #Parameters
    mu_max = 0.41
    Pob = 50.0
    Pma= 217.0
    Pmb = 108.0
    Pme = 127.0
    Ks = 0.5
    Ki = 200.0
    Si = 80.0
    Kmp = 0.5

    lmbda = 21.05
    Qpmax = 2.613
    Yps = 0.495
    beta = 0.0366

    delta = 0.8241
    alpha = 8.77
    b = 1.415
    a = 0.3142


    #Auxiliary Variables
    f_mu = 0.5*(1 - ((np.exp(lmbda*Z - delta) - np.exp(-lmbda*Z + delta))/(np.exp(lmbda*Z - delta) + np.exp(-lmbda*Z + delta))))

    P_ = (P- Pob)/(Pmb - Pob)
    if P <= Pob:
        P_ = 0
    elif P >= Pmb:
        P_ = 1

    S_ = S - Si
    if S <= Si:
        S_ = 0

    mu_ = (mu_max*S*(1 - (P/Pma)**a)*(1 - P_**b))/(Ks + S + (S*S_)/(Ki-Si))

    mu = f_mu*mu_

    Qp = Qpmax*(S/(Kmp + S))*(1-(P/Pme)**alpha)

    #define equations
    dXdt = (mu-D)*X

    dSdt = -(1/Yps)*Qp*X + D*Sf - D*S

    dPdt = Qp*X - D*P

    dZdt = beta*(W-Z)

    dWdt = beta*(Qp*X - D*P - W)

    return [dXdt, dSdt, dPdt, dZdt, dWdt]

def obj_calculator(sol, u, control_plan, kc, r, plot=False):
    #Extract t and states
    t_vals = sol.t
    x = sol.y

    #State_assignment
    state_assignment = {
        "X":0,
        "S":1,
        "P":2
    }

    if control_plan[1] == "Sf":
        Sf = np.zeros(len(t_vals))
        for i in range(len(t_vals)):
            Sf_ = kc*x[state_assignment[control_plan[0]]][i] + r
            if Sf_ < 0.0:
                Sf[i] = 0.0
            elif Sf_ > 200.0:
                Sf[i] = 200.0
            else:
                Sf[i] = Sf_
        D = u[0]*np.ones(len(t_vals))
    elif control_plan[1] == "D":
        D = np.zeros(len(t_vals))
        for i in range(len(t_vals)):
            D_ = kc*x[state_assignment[control_plan[0]]][i] + r
            if D_ < 0.0:
                D[i] = 0.0
            elif D_ > 0.06:
                D[i] = 0.06
            else:
                D[i] = D_
        Sf = u[1]*np.ones(len(t_vals))

    if plot == True:
        return t_vals, x, D, Sf


"""
Run controllers
"""
nominal_u = [0.06,200]
x0 = [10.0,100.0,0.0,0.0,0.0]
control_plan1 = ["P","D"]
kc1 = -2.5e-4
r1 = 0.04
sol_s_d = solve_ivp(model,[0,2500],x0,args=(nominal_u,control_plan1,kc1,r1),method="BDF")
tvals_s_d, x_s_d, D_s_d, Sf_s_d = obj_calculator(sol_s_d,nominal_u,control_plan1,kc1,r1,plot=True)

#Plot S_d
fig1,ax1, = plt.subplots(1,3, figsize=(8,3), num=1)
ax1[0].plot(tvals_s_d,10*x_s_d[0],label=r"10X")
ax1[0].plot(tvals_s_d,x_s_d[1],label=r"S")
ax1[0].plot(tvals_s_d,x_s_d[2],label=r"P")
ax1[0].plot(tvals_s_d,x_s_d[3],label=r"W")
ax1[0].plot(tvals_s_d,x_s_d[4],label=r"Z")
#Inputs 
ax1[1].plot(tvals_s_d,D_s_d)
ax1[2].plot(tvals_s_d,Sf_s_d)
#Formatting
ax1[0].legend(labelspacing=0.2,ncol=2, columnspacing=0.5, handlelength=0.5)
ax1[0].set_xlabel(r"Time (hr)")
ax1[1].set_xlabel(r"Time (hr)")
ax1[2].set_xlabel(r"Time (hr)")
ax1[0].set_ylabel(r"States")
ax1[1].set_ylabel(r"$D$ (1/hr)")
ax1[2].set_ylabel(r"$S_{f}$ (g/L)")
fig1.tight_layout()
plt.savefig("P_D_feedback_4.png",dpi=300)

plt.show()
