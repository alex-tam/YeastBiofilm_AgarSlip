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

"Data structure for dimensional parameters"
@with_kw struct Dimensional
    T::Float64 = 30420 # [min] Experiment duration
    Xc::Float64 = 2.25 # [mm] Initial biofilm half-width
    Xp::Float64 = 50.0 # [mm] Petri dish half-width
    ψn::Float64 # [mm^2/g/min] Cell profileration rate
    ψd::Float64 # [/min] Cell death rate
    G::Float64 = 5e-5 # [g/mm^2] Initial nutrient concentration
    D0::Float64 = 4.04e-2 # [mm^2/min] Nutrient diffusivity
    η::Float64 # [/min] Nutrient consumption rate
    Q::Float64 = 2.92e-3 # [mm/min] Nutrient mass transfer coefficient
    λ::Float64 # [g/mm^2/min] Slip coefficient
    μ::Float64 = 0.0534 # [g/mm/min] Biofilm (water) dynamic viscosity
    γ::Float64 = 0.0 # [g/mm/min^2] Surface tension coefficient
    ε::Float64 = 0.1 # [-] Thin-film parameter (aspect ratio)
    a::Float64 = 0.6 # [-] Agar density (%)
end

"Data structure for dimensionless parameters"
@with_kw struct Params
    ε::Float64 = 1e-6 # [-] Small parameter for Newton's method
    Nξ::Int = 101 # [-] Number of grid points
    Nτ::Int = 211 # [-] Number of time points
    plot_interval::Int = 210 # [-] Time points between output files
    S0::Float64 = 1 # [-] Initial contact line position
    H0::Float64 # [-] Initial biofilm height
    T::Float64 # [-] End time
    L::Float64 # [-] Domain width
    Ψd::Float64 # [-] Cell death rate
    Ψm::Float64 # [-] ECM production rate
    D::Float64 # [-] Nutrient diffusion coefficient
    Pe::Float64 # [-] Péclet number
    Υ::Float64 # [-] Nutrient consumption rate
    Qs::Float64 # [-] Nutrient loss rate
    Qb::Float64 # [-] Nutrient uptake rate
    λ::Float64 # [-] Slip coefficient
    γ::Float64 # [-] Surface tension coefficient
end

"Apply nondimensionalisation to obtain dimensionless parameters"
function nondimensionalise(dp, Ds, Db, ψm)
    # Apply nondimensionalisation
    nond_T = dp.ψn*dp.G*dp.T
    nond_L = dp.Xp/dp.Xc
    nond_Ψd = dp.ψd/(dp.ψn*dp.G)
    nond_Ψm = ψm/dp.ψn
    nond_D = Ds/(dp.G*dp.ψn*dp.Xc^2)
    nond_Pe = dp.ψn*dp.Xc^2*dp.G/Db
    nond_Υ = dp.η*dp.Xc/Db
    nond_Qs = dp.Q*dp.Xc/(dp.ε*Ds)
    nond_Qb = dp.Q*dp.Xc/(dp.ε*Db)
    nond_λ = dp.λ*dp.Xc/(dp.ε*dp.μ)
    nond_γ = dp.ε*dp.γ/(dp.ψn*dp.G*dp.Xc*dp.μ)
    par = Params(H0 = dp.ε, T = nond_T, L = nond_L, Ψd = nond_Ψd, Ψm = nond_Ψm, D = nond_D, Pe = nond_Pe, Υ = nond_Υ, Qs = nond_Qs, Qb = nond_Qb, λ = nond_λ, γ = nond_γ) # Initialise data structure for dimensionless parameters
    return par
end

"Main function for computing one solution"
function main()
    dp = Dimensional(ψn = 5.1, ψd = 2.57e-5, λ = 6.12e-5, η = 5.0e-4) # Load dimensional parameter data structure
    Ds = (1-0.023*dp.a)*dp.D0 # Nutrient diffusivity (agar)
    Db = 0.24*dp.D0 # Nutrient diffusivity (biofilm)
    ψm = dp.ψn/9 # ECM production rate
    par = nondimensionalise(dp, Ds, Db, ψm)
    println(par)
    # thinfilm_slip(par)
end

@time main()