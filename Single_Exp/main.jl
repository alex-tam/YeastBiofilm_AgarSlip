# Solve thin-film biofilm model with slip using front-fixing method
# Alex Tam, 12/01/2024

# Load packages
using Parameters
using LinearAlgebra
using Interpolations
using Printf
using Plots
using Measures
using LaTeXStrings

# Load external files
include("thinfilm_slip.jl")
include("solve_u.jl")
include("solve_h.jl")
include("solve_phi.jl")
include("solve_gs.jl")
include("solve_gb_ecm_advection.jl")
include("solve_S.jl")

"Data structure for dimensional parameters"
@with_kw struct Dimensional
    T::Float64 = 30420 # [min] Experiment duration
    Xc::Float64 = 2.25 # [mm] Initial biofilm half-width
    Xp::Float64 = 50.0 # [mm] Petri dish half-width
    ψn::Float64 # [mm^2/g/min] Cell profileration rate
    ψd::Float64 # [/min] Cell death rate
    G::Float64 = 5e-5 # [g/mm^2] Initial nutrient concentration
    D0::Float64 = 4.04e-2 # [mm^2/min] Glucose diffusivity in water
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

"Data structure for experiments"
struct ExpData
    a::Float64 # [-] Agar weight percentage
    t::Vector{Float64} # [min] Time for biofilm width collection
    w::Vector{Float64} # [mm] Vector of biofilm widths
    ϕ::Vector{Float64} # [-] Cell viability
    ar::Float64 # [-] Aspect ratio
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
    nond_Υ = dp.η*dp.Xc^2/Db
    nond_Qs = dp.Q*dp.Xc/(dp.ε*Ds)
    nond_Qb = dp.Q*dp.Xc/(dp.ε*Db)
    nond_λ = dp.λ*dp.Xc/(dp.ε*dp.μ)
    nond_γ = dp.ε*dp.γ/(dp.ψn*dp.G*dp.Xc*dp.μ)
    par = Params(H0 = dp.ε, T = nond_T, L = nond_L, Ψd = nond_Ψd, Ψm = nond_Ψm, D = nond_D, Pe = nond_Pe, Υ = nond_Υ, Qs = nond_Qs, Qb = nond_Qb, λ = nond_λ, γ = nond_γ) # Initialise data structure for dimensionless parameters
    return par
end
"Obtain experimental data for give agar density"
function get_exp(a)
    if a == 0.6 # 0.6% agar
        t = 24*60*[2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 21.0] # Experimental time
        w = 0.5*10*[0.5744081799793557, 0.9604482342331299, 1.2703540247559455, 1.5570392434668194, 1.802266346412321, 2.074337163491557, 2.2403945152741334, 2.7896823170919176]
        ϕ = 0.9*[85.85, 87.35, 91.75]/100 # Day 14 viability, converted to cell volume fraction
        ar = 2*0.019108 # [-] Experimental aspect ratio
    elseif a == 0.8 # 0.8% agar
        t = 24*60*[2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 21.0] # Experimental time
        w = 0.5*10*[0.6698935187956951, 1.026876275539356,1.3066850930981668, 1.5770597356826002, 1.825431441978235, 2.0917898707728235, 2.266949659727405, 2.8451534672347805]
        ϕ = 0.9*[85.85, 87.35, 91.75]/100 # Day 14 viability, converted to cell volume fraction
        ar = 2*0.019108*(1.72/1.36) # [-] Experimental aspect ratio
    elseif a == 1.2 # 1.2% agar
        t = 24*60*[2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 21.0] # Experimental time
        w = 0.5*10*[0.6164509723847065,0.9135083950335632,1.1340801867844186,1.3580911761585213,1.5651093345126472,1.7970215270038339,1.9523274477259278,2.475594028165591]
        ϕ = 0.9*[85.85, 87.35, 91.75]/100 # Day 14 viability, converted to cell volume fraction
        ar = 2*0.019108*(2.61/1.36) # [-] Experimental aspect ratio
    elseif a == 2.0 # 2.0% agar
        t = 24*60*[2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 21.0] # Experimental time
        w = 0.5*10*[0.5108086974571971,0.7167935742476953,0.869374371151493,1.0288330599608337,1.1677623597193019,1.367712440366207,1.4881608633666543,1.95958153860255]
        ϕ = 0.9*[85.85, 87.35, 91.75]/100 # Day 14 viability, converted to cell volume fraction
        ar = 2*0.019108*(2.92/1.36) # [-] Experimental aspect ratio
    else
        @printf("Invalid agar density.")
    end
    return t, w, ϕ, ar
end

"Main function for computing one solution"
function main()
    dp = Dimensional(ψn = 8.4, ψd = 1.6018e-5, λ = 0.00075, η = 0.9e-3) # Load dimensional parameter data structure
    Ds = (1-0.023*dp.a)*dp.D0 # Nutrient diffusivity (agar)
    Db = 0.24*dp.D0 # Nutrient diffusivity (biofilm)
    ψm = dp.ψn/9 # ECM production rate
    par = nondimensionalise(dp, Ds, Db, ψm)
    println(par)
    a = 0.6 
    t, w, ϕ, ar = get_exp(a)
    ex = ExpData(a, t, w, ϕ, ar) # Load experimental data
    thinfilm_slip(par, dp, ex, true)
end

@time main()