from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

#Plot formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size="14")
plt.rc('legend', fontsize=12)
plt.rc("axes", titlesize=12)


def model(t,x,u):
    X, S, P, Z, W = x

    #Inputs
    D,Sf = u

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

sol1 = solve_ivp(model, [0,800], [10.0, 100.0,0.0, 0.0,0.0],args=([0.04,200],),method="LSODA")
sol2 = solve_ivp(model, [0,800], [10.0, 100.0,0.0, 0.0,0.0],args=([0.12,200],),method="LSODA")
sol3 = solve_ivp(model, [0,800], [10.0, 100.0,0.0, 0.0,0.0],args=([0.06,50],),method="LSODA")

fig1,ax1 = plt.subplots(1,3,figsize=(8,3),num=1)
#Sol1
ax1[0].plot(sol1.t, 10*sol1.y[0],label=r"$10X$")
ax1[0].plot(sol1.t, sol1.y[1],label=r"$S$")
ax1[0].plot(sol1.t, sol1.y[2],label=r"$P$")
ax1[0].plot(sol1.t, sol1.y[3],label=r"$Z$")
ax1[0].plot(sol1.t, sol1.y[4],label=r"$W$")
# ax1[0].legend(ncol=2,labelspacing=0.2,handlelength=0.5)
ax1[0].set_xlabel(r"Time (hr)")
ax1[0].set_ylabel(r"States")
ax1[0].set_title(r"$D = 0.06, S_{f} = 200$")
#Sol2
ax1[1].plot(sol2.t, 10*sol2.y[0],label=r"$10X$")
ax1[1].plot(sol2.t, sol2.y[1],label=r"$S$")
ax1[1].plot(sol2.t, sol2.y[2],label=r"$P$")
ax1[1].plot(sol2.t, sol2.y[3],label=r"$Z$")
ax1[1].plot(sol2.t, sol2.y[4],label=r"$W$")
ax1[1].legend(ncol=2,labelspacing=0.2,handlelength=0.5)
ax1[1].set_xlabel(r"Time (hr)")
ax1[1].set_ylabel(r"States")
ax1[1].set_title(r"$D = 0.12, S_{f} = 200$")
#Sol3
ax1[2].plot(sol3.t, 10*sol3.y[0],label=r"$10X$")
ax1[2].plot(sol3.t, sol3.y[1],label=r"$S$")
ax1[2].plot(sol3.t, sol3.y[2],label=r"$P$")
ax1[2].plot(sol3.t, sol3.y[3],label=r"$Z$")
ax1[2].plot(sol3.t, sol3.y[4],label=r"$W$")
# ax1[2].legend(ncol=2,labelspacing=0.2,handlelength=0.5)
ax1[2].set_xlabel(r"Time (hr)")
ax1[2].set_ylabel(r"States")
ax1[2].set_title(r"$D = 0.06, S_{f} = 100$")
fig1.tight_layout()

# plt.savefig("ethanol_openloop.png",dpi=300)



plt.show()