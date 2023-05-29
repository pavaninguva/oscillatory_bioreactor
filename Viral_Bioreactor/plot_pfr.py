import pandas as pd
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

#Read simulation output
df = pd.read_csv("./pfr.csv")

#Plot
fig1, ax1 = plt.subplots(1,1,figsize=(5,4),num=1)
ax1.semilogy(df["tau"][:-15], df["T"][:-15],label=r"$T$")
ax1.semilogy(df["tau"][:-15], df["Id"][:-15],label=r"$I_{d}$")
ax1.semilogy(df["tau"][:-15], df["Is"][:-15],label=r"$I_{s}$")
ax1.semilogy(df["tau"][:-15], df["Ic"][:-15],label=r"$I_{c}$")
ax1.semilogy(df["tau"][:-15], df["Vs"][:-15],label=r"$V_{s}$")
ax1.semilogy(df["tau"][:-15], df["Vd"][:-15],label=r"$V_{d}$")
ax1.set_xlabel(r"$\tau$ (hr)")
ax1.set_ylabel(r"Concentration")
ax1.legend(ncol=3,handlelength=0.5)
plt.tight_layout()

plt.savefig("frensing_pfr.png",dpi=300)

plt.show()