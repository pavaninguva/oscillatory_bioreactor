from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

#Plot formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size="14")


def model (t,x,u,params):
    #Unpack x, u, params
    T = x[0]
    Id = x[1]
    Is = x[2]
    Ic = x[3]
    Vs = x[4]
    Vd  = x[5]

    Tin = u[0]
    D = u[1]

    mu = params[0]
    k1 = params[1]
    k2 = params[2]
    k3 = params[3]
    k33 = params[4]
    k4 = params[5]
    f = params[6]

    #Define model
    dTdt = mu*T - k1*(Vs+Vd)*T + (D)*(Tin-T)
    dIddt = k1*Vd*T - (k1*Vs - mu)*Id - (D)*Id 
    dIsdt = k1*Vs*T -(k1*Vd + k2)*Is -(D)*Is 
    dIcdt = k1*(Vs*Id + Vd*Is) -k2*Ic - (D)*Ic 
    dVsdt = k3*Is - (k1*(T+Id+Is+Ic) + k4 + (D))*Vs 
    dVddt = k33*Ic + f*k3*Is - (k1*(T+Id+Is+Ic) + k4 + (D))*Vd 

    return [dTdt,dIddt,dIsdt,dIcdt,dVsdt,dVddt]


nominal_params = [0.027,2.12e-9,7.13e-3,168,168,0.035,1e-3]
nominal_u1 = [3e6,0.0396]
nominal_u2 = [3e6,1.05]
nominal_u3 = [3e6,0.25]
x0 = [5e6,0.0,0.0,0.0,1.25e5,0.0]


sol1 = solve_ivp(model,[0,400],x0,args=(nominal_u1,nominal_params))
sol2 = solve_ivp(model,[0,400],x0,args=(nominal_u2,nominal_params))
sol3 = solve_ivp(model,[0,400],x0,args=(nominal_u3,nominal_params))


#Plot
fig1, ax1 = plt.subplots(1,1,figsize=(4,3),num=1)
ax1.semilogy(sol1.t,sol1.y[0],label=r"$T$")
ax1.plot(sol1.t,sol1.y[1],label=r"$I_{d}$")
ax1.plot(sol1.t,sol1.y[2],label=r"$I_{s}$")
ax1.plot(sol1.t,sol1.y[3],label=r"$I_{c}$")
ax1.plot(sol1.t,sol1.y[4],label=r"$V_{s}$")
ax1.plot(sol1.t,sol1.y[5],label=r"$V_{d}$")
ax1.set_xlabel(r"Time (hr)")
ax1.set_ylabel(r"Concentration")
ax1.legend(ncol=2,columnspacing=0.7)
fig1.tight_layout()
plt.savefig("nominal.png",dpi=300)

fig2, ax2 = plt.subplots(1,1,figsize=(4,3),num=2)
ax2.semilogy(sol2.t,sol2.y[0],label=r"$T$")
ax2.plot(sol2.t,sol2.y[1],label=r"$I_{d}$")
ax2.plot(sol2.t,sol2.y[2],label=r"$I_{s}$")
ax2.plot(sol2.t,sol2.y[3],label=r"$I_{c}$")
ax2.plot(sol2.t,sol2.y[4],label=r"$V_{s}$")
ax2.plot(sol2.t,sol2.y[5],label=r"$V_{d}$")
ax2.set_xlabel(r"Time (hr)")
ax2.set_ylabel(r"Concentration")
ax2.legend(ncol=2,columnspacing=0.7)
fig2.tight_layout()
plt.savefig("washout.png",dpi=300)

fig3, ax3 = plt.subplots(1,1,figsize=(4,3),num=3)
ax3.semilogy(sol3.t,sol3.y[0],label=r"$T$")
ax3.plot(sol3.t,sol3.y[1],label=r"$I_{d}$")
ax3.plot(sol3.t,sol3.y[2],label=r"$I_{s}$")
ax3.plot(sol3.t,sol3.y[3],label=r"$I_{c}$")
ax3.plot(sol3.t,sol3.y[4],label=r"$V_{s}$")
ax3.plot(sol3.t,sol3.y[5],label=r"$V_{d}$")
ax3.set_xlabel(r"Time (hr)")
ax3.set_ylabel(r"Concentration")
ax3.legend(ncol=2,columnspacing=0.7)
fig3.tight_layout()
plt.savefig("stable.png",dpi=300)

plt.show()