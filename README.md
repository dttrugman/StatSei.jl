# StatSei.jl
Statistical Seismology Tools in Julia


This repository hosts various Julia scripts for performing statistical seismology tasks. It is currently a work in progress and under development, so please have patience.

---

The (unregistered) package can be installed using the Julia Pkg manager:

` pkg> add https://github.com/dttrugman/StatSei.jl`

The example notebooks use additional functionality from the following packages, so if you would like to run them please also add the following:

` pkg> add DataFrames, GLM, CSV, StatsBase, PyPlot`

Note that at present writing, the plotting functionality involves `PyPlot` which while functional is not ideal for a Julia package. Pull requests to help fix this are most welcome.

---

The `notebooks/` directory contains several example notebooks doing various calculations:
