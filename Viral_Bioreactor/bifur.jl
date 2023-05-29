using Revise, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

# define the sup norm
norminf(x) = norm(x, Inf)

# function to record information from a solution
recordFromSolution(x, p) = (T = x[1], Id = x[2], Is = x[3], 
							Ic = x[4], Vs=x[5], Vd=x[6])

function pp2!(dz, z, p, t)
	@unpack Tin, D, mu, k1, k2, k3, k33, k4, f = p
	T, Id, Is, Ic, Vs, Vd = z
	dz[1] = mu*T - k1*(Vs+Vd)*T + (D)*(Tin-T)
	dz[2] =	k1*Vd*T - (k1*Vs - mu)*Id - (D)*Id
    dz[3] = k1*Vs*T -(k1*Vd + k2)*Is -(D)*Is
    dz[4] = k1*(Vs*Id + Vd*Is) -k2*Ic - (D)*Ic
    dz[5] = k3*Is - (k1*(T+Id+Is+Ic) + k4 + (D))*Vs
    dz[6] = k33*Ic + f*k3*Is - (k1*(T+Id+Is+Ic) + k4 + (D))*Vd
	dz
end

pp2(z, p) = pp2!(similar(z), z, p, 0)

# parameters of the model
par_pp2 = (Tin = 3e6, D = 0.0396, mu = 0.027, k1 = 2.12e-9,
			k2 = 7.13e-3, k3 = 168, k33 = 168, k4 = 0.035, f = 1e-3)

#initial conditions
z0 = [5e6,1e6,1e6,1e6,1e6,1e7]

# bifurcation problem
prob = BifurcationProblem(pp2, z0, par_pp2,
	# specify the continuation parameter
	(@lens _.D), recordFromSolution = recordFromSolution)

# continuation options
opts_br = ContinuationPar(pMin = 0.03, pMax = 1.2, dsmax = 0.01,
	# options to detect bifurcations
	detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 50,
	# number of eigenvalues
	nev = 6,
	# maximum number of continuation steps
	maxSteps = 1000,)

diagram = bifurcationdiagram(prob, PALC(),
	3,
	(args...) -> setproperties(opts_br; ds = -0.001, dsmax = 0.01, nInversion = 8, detectBifurcation = 3);
	# δp = -0.01,
	verbosity = 0, plot = true)

scene = plot(diagram; code = (), title="$(size(diagram)) branches", legend = false)


