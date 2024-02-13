module StatSei

# external packages
using PyPlot
using StatsBase
using DataFrames
using Dates
using CSV
using Printf

# and modules within specific packages
using Distances: pairwise, Haversine, Euclidean
using Proj: Transformation
using Interpolations: LinearInterpolation, linear_interpolation
using GaussianMixtures: GMM
using Distributions: Normal

# exports from catalogs
include("catalogs.jl")
export load_catalog_scsn, load_catalog_hs2011, load_catalog_hs2022
export load_catalog_ncsn, load_catalog_nsl
export load_catalog_nvreloc, load_catalog_merged

# exports from fmd
include("fmd.jl")
export fmd, mc_maxc, bval_maxl, bval_bpos

# exports from neighbors
include("neighbors.jl")
export calculate_erad, calculate_l10Rad, calculate_l10LWC
export find_neighbors_llt, find_neighbors_xyt, find_neighbors_xyzt
export get_thresh_shuffle_llt, get_thresh_shuffle_xyt, get_thresh_shuffle_xyzt
export fit_gmm, assign_families, plot_neighbors

# exports from pairwise
include("pairwise.jl")
export find_neighbors_pllt, find_neighbors_pxyt, find_neighbors_pxyzt

# export from projection
include("projection.jl")
export lonlat2xypos, xypos2latlon, map_distance, xydist, setup_projection

end # module StatSei
