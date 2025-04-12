# Parameter estimation for thin-film slip model
# Alex Tam, 17/02/2025

# Load packages
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

# Load external files
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
    Xc::Float64 = 1.5 # [mm] Initial biofilm half-width
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
    ε::Float64 = 0.06 # [-] Thin-film parameter (aspect ratio)
    H0::Float64 = 0.1 # [-] Initial biofilm height
end

"Data structure for dimensionless parameters"
@with_kw struct Params
    ε::Float64 = 1e-6 # [-] Small parameter for Newton's method
    Nξ::Int = 101 # [-] Number of grid points
    Nτ::Int = 841 # [-] Number of time points
    plot_interval::Int = 210 # [-] Time points between output files
    S0::Float64 = 1 # [-] Initial contact line position
    T::Float64 # [-] End time
    L::Float64 # [-] Domain width
    Ψd::Float64 # [-] Cell death rate
    D::Float64 # [-] Nutrient diffusion coefficient
    Pe::Float64 # [-] Péclet number
    Υ::Float64 # [-] Nutrient consumption rate
    Qs::Float64 # [-] Nutrient loss rate
    Qb::Float64 # [-] Nutrient uptake rate
    λ::Float64 # [-] Slip coefficient
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
function nondimensionalise(dp, Ds, Db)
    # Apply nondimensionalisation
    nond_T = dp.ψn*dp.G*dp.T
    nond_L = dp.Xp/dp.Xc
    nond_Ψd = dp.ψd/(dp.ψn*dp.G)
    nond_D = Ds/(dp.G*dp.ψn*dp.Xc^2)
    nond_Pe = dp.ψn*dp.Xc^2*dp.G/Db
    nond_Υ = dp.η*dp.Xc^2/Db
    nond_Qs = dp.Q*dp.Xc/(dp.ε*Ds)
    nond_Qb = dp.Q*dp.Xc/(dp.ε*Db)
    nond_λ = dp.λ*dp.Xc/(dp.ε*dp.μ)
    par = Params(T = nond_T, L = nond_L, Ψd = nond_Ψd, D = nond_D, Pe = nond_Pe, Υ = nond_Υ, Qs = nond_Qs, Qb = nond_Qb, λ = nond_λ) # Initialise data structure for dimensionless parameters
    return par
end

"Alternative method for fine-grid solutions"
function nondimensionalise(dp, Ds, Db, ξ, τ, pl_int)
    # Apply nondimensionalisation
    nond_T = dp.ψn*dp.G*dp.T
    nond_L = dp.Xp/dp.Xc
    nond_Ψd = dp.ψd/(dp.ψn*dp.G)
    nond_D = Ds/(dp.G*dp.ψn*dp.Xc^2)
    nond_Pe = dp.ψn*dp.Xc^2*dp.G/Db
    nond_Υ = dp.η*dp.Xc/Db
    nond_Qs = dp.Q*dp.Xc/(dp.ε*Ds)
    nond_Qb = dp.Q*dp.Xc/(dp.ε*Db)
    nond_λ = dp.λ*dp.Xc/(dp.ε*dp.μ)
    par = Params(Nξ = ξ, Nτ = τ, plot_interval = pl_int, T = nond_T, L = nond_L, Ψd = nond_Ψd, D = nond_D, Pe = nond_Pe, Υ = nond_Υ, Qs = nond_Qs, Qb = nond_Qb, λ = nond_λ) # Initialise data structure for dimensionless parameters
    return par
end

"Compute the objective function: distance between model and parameters for given solution"
function objective(p::Vector{T}, ex::ExpData) where{T}
    ##### Obtain dimensionless parameters
    dp = Dimensional(ψn = p[1], ψd = p[2], η = p[3], λ = p[4]) # Load dimensional parameter data structure
    Ds = (1-0.023*ex.a)*dp.D0 # Nutrient diffusivity (agar)
    Db = 0.24*dp.D0 # Nutrient diffusivity (biofilm)
    par = nondimensionalise(dp, Ds, Db)
    ##### Compute model solution
    dist = thinfilm_slip(par, dp, ex, false)
    return dist
end

"Function for both the objective and gradient"
function fg!(F, G, p::Vector{T}, ex::ExpData) where{T}
    δ::Float64 = 1e-6
    ##### Compute objective function
    dp = Dimensional(ψn = p[1], ψd = p[2], η = p[3], λ = p[4]) # Load dimensional parameter data structure
    Ds = (1-0.023*ex.a)*dp.D0 # Nutrient diffusivity (agar)
    Db = 0.24*dp.D0 # Nutrient diffusivity (biofilm)
    par = nondimensionalise(dp, Ds, Db)
    dist = thinfilm_slip(par, dp, ex, false)
    ##### Compute gradient
    if G !== nothing
        # Gradient in ψn direction
        dp1 = Dimensional(ψn = p[1] + δ, ψd = p[2], η = p[3], λ = p[4]) # Load dimensional parameter data structure
        par1 = nondimensionalise(dp1, Ds, Db)
        p1_dist = thinfilm_slip(par1, dp1, ex, false)
        G[1] = (p1_dist-dist)/(δ)
        # Gradient in ψd direction
        dp2 = Dimensional(ψn = p[1], ψd = p[2] + δ, η = p[3], λ = p[4]) # Load dimensional parameter data structure
        par2 = nondimensionalise(dp2, Ds, Db)
        p2_dist = thinfilm_slip(par2, dp2, ex, false)
        G[2] = (p2_dist-dist)/(δ)
        # Gradient in η direction
        dp3 = Dimensional(ψn = p[1], ψd = p[2], η = p[3] + δ, λ = p[4]) # Load dimensional parameter data structure
        par3 = nondimensionalise(dp3, Ds, Db)
        p3_dist = thinfilm_slip(par3, dp3, ex, false)
        G[3] = (p3_dist-dist)/(δ)
        # Gradient in λ direction
        dp4 = Dimensional(ψn = p[1], ψd = p[2], η = p[3], λ = p[4] + δ) # Load dimensional parameter data structure
        par4 = nondimensionalise(dp4, Ds, Db)
        p4_dist = thinfilm_slip(par4, dp4, ex, false)
        G[4] = (p4_dist-dist)/(δ)
    end
    if F !== nothing
        return dist
    end
end

"Main function for parameter estimation"
function main()
    A = [0.6, 0.8, 1.2, 2.0] # Agar density
    Vp = Vector{Vector{Float64}}()
    Vf = Vector{Float64}()
    for a in A
        ##### Experimental data
        t, w, ϕ, ar = get_exp(a)
        ex = ExpData(a, t, w, ϕ, ar) # Load experimental data
        ##### Optimise parameter estimates for given experimental data
        # Global black box optimisation
        res = bboptimize(x -> objective(x, ex); SearchRange = [(5.0, 25.0), (0.0, 5e-4), (0.0, 1e-2), (0.0, 1e-2)], 
            MaxFuncEvals = 500, Method = :random_search)
        p0 = best_candidate(res)
        lower = [0.0, 0.0, 0.0, 0.0]
        upper = [Inf, Inf, Inf, Inf]
        # Gradient-based box-constrained optimisation
        inner_optimizer = Optim.LBFGS(m = 20, 
            linesearch = LineSearches.HagerZhang(),
            alphaguess = LineSearches.InitialHagerZhang(α0 = 1.0)) # Constructor for box-constrained solver
        @time result = Optim.optimize(Optim.only_fg!((F, G, x) -> fg!(F, G, x, ex)), 
            lower, upper, p0, 
            Fminbox(inner_optimizer), 
            Optim.Options(outer_iterations = 3, 
                x_abstol = 1e-6, f_abstol = 1e-6, g_abstol = 1e-6, 
                outer_x_abstol = 1e-6, outer_f_abstol = 1e-6, outer_g_tol = 1e-6, 
                show_trace = true))
        p = Optim.minimizer(result)
        f = Optim.minimum(result)
        push!(Vp, p)
        push!(Vf, f)
        writedlm("OptimalParameters.csv", Vp)
        writedlm("OptimalObjectives.csv", Vf)
        ##### Plot optimal solution
        dp = Dimensional(ψn = p[1], ψd = p[2], η = p[3], λ = p[4]) # Load dimensional parameter data structure
        Ds = (1-0.023*ex.a)*dp.D0 # Nutrient diffusivity (agar)
        Db = 0.24*dp.D0 # Nutrient diffusivity (biofilm)
        par = nondimensionalise(dp, Ds, Db)
        dist = thinfilm_slip(par, dp, ex, true)
    end
end

@time main()