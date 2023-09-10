# open contents of Project.toml (explicitly added packages) and Manifest.toml (recursive dependency packages)
import Pkg; Pkg.activate(".")

# # run the following if no Project.toml and Manifest.toml files are present in the current dir
# # will take a few minutes to run if first time creating full env
# Pkg.add("Conda") # for managing python packages made private to Julia
# Pkg.add("PyCall") # for calling python scripts inside Julia
# Pkg.add("JLD") # for loading and saving Julia variables using .jld file format
# Pkg.add("LaTeXStrings") # for using latex math fonts in plots, prepend L to strings, ie. L"5+5"
# Pkg.add("BenchmarkTools") # for performance tracking, https://github.com/JuliaCI/BenchmarkTools.jl
# Pkg.add("OrdinaryDiffEq") # for just the ODE numerical solvers in DifferentialEquations.jl
# Pkg.add("DiffEqCallbacks") # for saving callback
# Pkg.add("Documenter") # for documenting your functions with docstrings that support latex, see https://stackoverflow.com/questions/44120903/julia-docstrings-and-latex
# Pkg.add("Plots") # for plotting in Julia 
# Pkg.add("PyPlot") # matplotlib backend for Julia plotting

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
using Plots; pyplot(html_output_format=:png) # activates `pyplot()`` backend instead of `gr()` default
using Random

# show status of the pkg env
Pkg.status()

# need to set to empty string for python dir and recompile to force Julia installation 
# to use own mini Python distribution, more info at: https://github.com/JuliaPy/PyCall.jl
ENV["PYTHON"]=""
println("Conda.jl takes packages from $(Conda.ROOTENV)")
println("PyCall uses python at $(Conda.PYTHONDIR)")
println("More precisely PyCall runs using packages at $(PyCall.pyprogramname)")
Conda.list(Conda.ROOTENV)

# # install some packages into Julia's python env, only need to run once
# Conda.add("numpy") # scientific computing package
# Conda.add("scipy") # scientific computing package
# Conda.add("clawpack") # clawpack only works in linux / macOS, use WSL2 in Windows

# test julia
RNG = Random.seed!(MersenneTwister(1))
testtest = rand(RNG, 1) # should be reproducible random number at this point
formatted_label = "$(@sprintf("%.2f", testtest[1]))" # pretty formatting
println(formatted_label)
println(L"5+5")
date = "07-08-2018"
titlestring = latexstring("\\mathrm{My\\, \\alpha\\, including\\, a\\, variable\\, date:\\,}", date)
println(titlestring)
display(L"\textrm{My date is } %$(date)")

# need to run in Jupyer notebook to visualize
plot(
    sin,
    (0:0.01:5*pi),
    title = titlestring,
    xlabel = L"\mathrm{My\, x-label\, (mJ/cm^{3})}",
    ylabel = L"\mathrm{\alpha (m/s)}",
     label = latexstring("5x+4y")
)

# Usage of file I/O via JLD.jl
jld_fp = pwd() * "/TestTest.jld"
a = Dict{String, Array}("A"=>[1,23], "B"=>[5.5,2], "15"=>[0;0;0])

jldopen(jld_fp, "w") do file
    write(file, "a", a)  # alternatively: "@write file a"
end

c = jldopen(jld_fp, "r") do file
    read(file, "a")
end
display(c)
rm(jld_fp)

# test python using pycall
py"""
# packages for system vars
import sys
import os

print(os.system("pwd"))
print(sys.version)
"""

# test PyClaw, solve acoustics equations, more info at:
# http://www.clawpack.org/pyclaw/tutorial.html#pyclaw-tutorial
py"""
# import relevant packages
import numpy as np
from scipy.stats import norm
from clawpack import pyclaw
from clawpack import riemann
from math import sqrt
from numpy import exp

# solver
solver = pyclaw.ClawSolver1D(riemann.acoustics_1D)
solver.bc_lower[0] = pyclaw.BC.wall
solver.bc_upper[0] = pyclaw.BC.extrap

# domain
domain = pyclaw.Domain([-1.0], [1.0], [200])
solution = pyclaw.Solution(solver.num_eqn, domain)

# initial condition
state = solution.state
xc = state.grid.p_centers[0]
state.q[0,:] = exp(-100 * (xc-0.75)**2)
state.q[1,:] = 0.0

# parameters
rho = 1.0
bulk = 1.0
state.problem_data['rho'] = rho
state.problem_data['bulk'] = bulk
state.problem_data['zz'] = sqrt(rho*bulk)
state.problem_data['cc'] = sqrt(bulk/rho)

# Controller
controller = pyclaw.Controller()
controller.solution = solution
controller.solver = solver
controller.tfinal = 1.0
controller.keep_copy = True

# run
status = controller.run()
"""

# Just use Julia pyplot() backend (==matplotlib() wrapper) to plot the desired results
# need to run in Jupyer notebook to visualize
centers = py"controller.centers"
results = py"controller.frames"
plot(
    centers,
    results[end].q[1,:] 
)