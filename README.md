# StatSei.jl
Statistical Seismology Tools in Julia


This repository hosts various Julia scripts for performing statistical seismology tasks. It is currently a work in progress and under development, so please have patience.

This compilation of codes are dedicated to Ilya Zaliapin, a dear friend, a mentor, and an outstanding scientist whose ideas will resonate through the statistical seismology community for years to come.

---

The (unregistered) package can be installed using the Julia Pkg manager:

` pkg> add https://github.com/dttrugman/StatSei.jl`


Note that at present writing, the plotting functionality involves `PyPlot` which while functional is not ideal for a Julia package. Pull requests to help fix this are most welcome.

---

The `notebooks/` directory contains several example notebooks doing various calculations. These notebooks use additional functionality from the external packages, so if you would like to run them please also add the following:

` pkg> add DataFrames, GLM, CSV, StatsBase, PyPlot`

The examples are as follows:

- `fmd_bval_example.ipynb`: example calculations involving frequency-magnitude distributions (FMDs) and b-value estimation.

- `nearest_neighbors_example.ipynb`: example implementation of the nearest neighbor seismicity clustering analysis proposed by Zaliapin et al. (2008) and Zaliapin and Ben-Zion (2013).

- `projection_fractal_example.ipynb`: examples of transformation of seismicity catalogs to cartesian coordinates in order to calculate the spatial fractal dimension.
