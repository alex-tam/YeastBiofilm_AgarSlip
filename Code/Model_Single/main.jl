# Solve thin-film slip model
# Alex Tam, 13/04/2025

# Import packages
using Revise
using Parameters
using LinearAlgebra
using Interpolations
using Printf
using DelimitedFiles
using Plots
using Measures
using LaTeXStrings
using BlackBoxOptim
using Optim
using LineSearches

# Load code files
include("get_exp.jl")
include("thinfilm_slip.jl")
include("solve_u.jl")
include("solve_h.jl")
include("solve_phi.jl")
include("solve_gs.jl")
include("solve_gb.jl")
include("solve_S.jl")

"Data structure for dimensional parameters"
@with_kw struct Dimensional
    T::Float64 = 30780 # [min] Experiment duration
    Hs::Float64 = 3.0 # [mm] Substratum depth
    Xs::Float64 = 50.0 # [mm] Substratum length scale (Petri dish half-width)
    Xb::Float64 = 2.0 # [mm] Biofilm length scale (Initial biofilm half-width)
    H0::Float64 = 0.006 # [mm] Initial biofilm height
    G::Float64 = 5e-5 # [g/mm^2] Initial nutrient concentration
    D0::Float64 = 4.04e-2 # [mm^2/min] Glucose diffusivity in water
end

"Data structure for dimensionless parameters"
@with_kw struct Params
    δ::Float64 = 1e-6 # [-] Small parameter for Newton's method
    Nξ::Int = 101 # [-] Number of grid points (biofilm)
    Nξo::Int = 101 # [-] Number of grid points (substratum)
    Nτ::Int = 401 # [-] Number of time points
    plot_interval::Int = 40 # [-] Time points between output files
    T::Float64 # [-] End time
    L::Float64 # [-] Domain width
    H0::Float64 # [-] Initial biofilm height
    S0::Float64 = 1.0 # [-] Initial contact line position
    ε::Float64 # [-] Substratum aspect ratio
    Ψn::Float64 = 0.30224739765610387 # [-] Biomass proliferation rate
    Ψd::Float64 = 0.008345663137850129 # [-] Biomass death rate
    D::Float64 # [-] Nutrient diffusion coefficient
    Q::Float64 = 3.272246675181909 # [-] Nutrient uptake rate
    Υ::Float64 = 5.010540755996383 # [-] Nutrient consumption rate
    λ::Float64 = 2.6707602643634694 # [-] Slip coefficient
end

"Data structure for experimental data"
struct ExpData
    a::Float64 # [-] Agar weight percentage
    t::Vector{Float64} # [-] Time for biofilm width collection
    w::Vector{Float64} # [-] Vector of biofilm widths
    ϕ::Vector{Float64} # [-] Cell viability
    ar::Float64 # [-] Aspect ratio
end

"Apply nondimensionalisation to obtain dimensionless parameters"
function nondimensionalise(dp, Ds, Db)
    # Apply nondimensionalisation
    non_T = Db/dp.Xb^2*dp.T
    non_L = dp.Xs/dp.Xb
    non_H0 = dp.H0*dp.Xs/(dp.Hs*dp.Xb)
    non_ε = dp.Hs/dp.Xs
    non_D = Ds/Db
    return Params(T = non_T, L = non_L, H0 = non_H0, ε = non_ε, D = non_D)  				
end

"Main function for solving thin-film model"
function main()
    ##### Establish plots
    gr() # Load GR plotting back-end
    default(legendfont = (16, "Computer Modern"), guidefont = (24, "Computer Modern"), tickfont = (20, "Computer Modern")) # Plot settings
    ##### Initialise parameters and experimenta data
    dp = Dimensional()
    a = 2.0
    Db = 0.24*dp.D0 # [mm^2/min] Nutrient diffusivity (biofilm)
    Ds = (1-0.023*a)*dp.D0
    par = nondimensionalise(dp, Ds, Db)
    t, w, ϕ, ar = get_exp(a, dp, Db)
    ex = ExpData(a, t, w, ϕ, ar) # Load experimental data
    ##### Solve model
    dist = thinfilm_slip(par, ex, true)
    return dist
end

@time main()