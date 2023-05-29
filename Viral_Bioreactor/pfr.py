from fipy import Grid1D, CellVariable, TransientTerm, ExplicitUpwindConvectionTerm, Viewer, ImplicitSourceTerm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

nx = 100
Lx = 20.0
mesh = Grid1D(nx=nx, Lx = Lx)

T = CellVariable(mesh=mesh, hasOld=True, value = 0.0,name="T")
Id = CellVariable(mesh=mesh, hasOld=True, value=0.0,name="I_d")
Is = CellVariable(mesh=mesh, hasOld=True, value=0.0,name="I_s")
Ic = CellVariable(mesh=mesh, hasOld=True, value=0.0,name="I_c")
Vs = CellVariable(mesh=mesh, hasOld=True, value= 0.0,name="V_s")
Vd = CellVariable(mesh=mesh, hasOld=True, value=0.0,name="V_d")

mu, k1, k2, k3, k33, k4, f = [0.027,2.12e-9,7.13e-3,168,168,0.035,1e-3]


#Boundary conditions
T0, Id0, Is0, Ic0, Vs0, Vd0 = [5e6,0.0,0.0,0.0,1.25e5,0.0]

#Inlet
T.constrain(T0, mesh.facesLeft)
Id.constrain(Id0, mesh.facesLeft)
Is.constrain(Is0, mesh.facesLeft)
Ic.constrain(Ic0, mesh.facesLeft)
Vs.constrain(Vs0, mesh.facesLeft)
Vd.constrain(Vd0, mesh.facesLeft)

#Outlet
T.faceGrad.constrain((0),mesh.facesRight)
Id.faceGrad.constrain((0),mesh.facesRight)
Is.faceGrad.constrain((0),mesh.facesRight)
Ic.faceGrad.constrain((0),mesh.facesRight)
Vs.faceGrad.constrain((0),mesh.facesRight)
Vd.faceGrad.constrain((0),mesh.facesRight)

#Model Definition

# 6m/hr 
convcoeff = (1.0,)

eq_T = (TransientTerm(var=T) + ExplicitUpwindConvectionTerm(coeff=convcoeff, var=T)
        # + k1*(Vs+Vd)*T - mu*T
        -ImplicitSourceTerm(coeff=mu,var=T) + ImplicitSourceTerm(coeff=k1*(Vs+Vd),var=T)
) == 0

eq_Id = (TransientTerm(var=Id) + ExplicitUpwindConvectionTerm(coeff=convcoeff, var=Id)
        #  + k1*Vd*T + (k1*Vd - mu)*Id
        - ImplicitSourceTerm(coeff = k1*(Vd),var=T) + ImplicitSourceTerm(coeff=(k1*Vd - mu),var=Id)
) == 0

eq_Is = (TransientTerm(var=Is) + ExplicitUpwindConvectionTerm(coeff=convcoeff, var=Is)
        #  -k1*Vs*T + (k1*Vd + k2)*Is
        - ImplicitSourceTerm(coeff = k1*(Vs),var=T) + ImplicitSourceTerm(coeff=(k1*Vd + k2),var=Is)
) == 0


eq_Ic = (TransientTerm(var=Ic) + ExplicitUpwindConvectionTerm(coeff=convcoeff,var=Ic)
        #  + k2*Ic - k1*Vs*Id - k1*Vd*Is
         + ImplicitSourceTerm(coeff=k2,var=Ic) - ImplicitSourceTerm(coeff=(k1*Vs),var=Id) - ImplicitSourceTerm(coeff=(k1*Vd),var=Is)
) == 0

eq_Vs = (TransientTerm(var=Vs) + ExplicitUpwindConvectionTerm(coeff=convcoeff,var=Vs)
        #  -k3*Is + (k1*(T+Id+Is+Ic) + k4)*Vs
         -ImplicitSourceTerm(coeff=k3,var=Is) + ImplicitSourceTerm(coeff = (k1*(T+Id+Is+Ic) + k4),var=Vs)
) == 0

eq_Vd = (TransientTerm(var=Vd) + ExplicitUpwindConvectionTerm(coeff=convcoeff,var=Vd)
        #  - k33*Ic + (k1*(T+Id+Is+Ic) + k4)*Vd - f*k3*Is
         -ImplicitSourceTerm(coeff=k33,var=Ic) + ImplicitSourceTerm(coeff = (k1*(T+Id+Is+Ic) + k4),var=Vd) - ImplicitSourceTerm(coeff=(f*k3),var=Is)
) == 0

#Couple Equations
eq_main = eq_T & eq_Id & eq_Is & eq_Ic & eq_Vs & eq_Vd

#Initialize viewer
# if __name__ == "__main__":
    # viewer = Viewer((T, Id, Is, Ic, Vs, Vd),ylog=True)


#Solve
CFL = 1.0
dt = CFL * (Lx/nx)/convcoeff[0]

t = 0.0
t_end = 20.0
counter = 0
t_vals = []


while t < t_end:
    t_vals.append(t)
    t += dt
    counter += 1
    T.updateOld()
    Is.updateOld()
    Ic.updateOld()
    Id.updateOld()
    Vs.updateOld()
    Vd.updateOld()

    res = 1e5
    while res > 1e-5:
        res = eq_main.sweep(dt=dt)
        print (t)

    if t > 17.0 - 1e-8:
        df = pd.DataFrame({
            "tau":mesh.x,
            "T":T,
            "Is":Is,
            "Ic":Ic,
            "Id":Id,
            "Vs":Vs,
            "Vd":Vd
        })
        break


df.to_csv("pfr.csv")
        
    
"""
Export as csv
"""
