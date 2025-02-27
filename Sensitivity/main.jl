# Solve thin-film biofilm model with slip using front-fixing method
# Alex Tam, 12/01/2024

# Load packages
using Parameters
using LinearAlgebra
using Interpolations
using Printf
using DelimitedFiles

# Load external files
include("thinfilm_slip.jl")
include("solve_u.jl")
include("solve_h.jl")
include("solve_phi.jl")
include("solve_gs.jl")
include("solve_gb.jl")
include("solve_S.jl")

"Data structure for numerical and model parameters"
@with_kw struct Params
    T::Float64 = 10.0 # [-] End time
    Nτ::Int = 10001 # [-] Number of time points
    plot_interval::Int = 1000 # [-] Time points between output files
    L::Float64 = 10.0 # [-] Domain width
    Nξ::Int = 501 # [-] Number of grid points
    S0::Float64 = 1 # [-] Initial contact line position
    H0::Float64 = 0.1 # [-] Initial biofilm height
    Ψd::Float64 = 0.0 # [-] Cell death rate
    Ψm::Float64 = 0.05 # [-] ECM production rate
    D::Float64 = 4.34 # [-] Nutrient diffusion coefficient
    Pe::Float64 = 0.953 # [-] Péclet number
    Υ::Float64 = 3.15 # [-] Nutrient consumption rate
    Qs::Float64 = 2.09 # [-] Nutrient loss rate
    Qb::Float64 = 8.65 # [-] Nutrient uptake rate
    λ::Float64 = 0.0 # [-] Slip coefficient
    γ::Float64 = 0.0 # [-] Surface tension coefficient
    ε::Float64 = 1e-6 # [-] Small parameter for Newton's method
end

"Main function for varying model parameters"
function main()
    ### Configure solutions ###
    parameters = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05] # Parameter values
    ### Compute solutions ###
    for p in parameters
        par = Params(Ψd = p)
        @time thinfilm_slip(par, p)
    end
end

@time main()