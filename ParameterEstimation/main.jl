# Parameter estimation for thin-film slip model
# Alex Tam, 13/04/2025

# Import packages
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
    Xb::Float64 = 1.5 # [mm] Biofilm length scale (Initial biofilm half-width)
    H0::Float64 = 0.006 # [mm] Initial biofilm height
    G::Float64 = 5e-5 # [g/mm^2] Initial nutrient concentration
    D0::Float64 = 4.04e-2 # [mm^2/min] Glucose diffusivity in water
    Q::Float64 = 2.92e-3 # [mm/min] Nutrient mass transfer coefficient
end

"Data structure for dimensionless parameters"
@with_kw struct Params
    δ::Float64 = 1e-6 # [-] Small parameter for Newton's method
    Nξ::Int = 101 # [-] Number of grid points (biofilm)
    Nξo::Int = 101 # [-] Number of grid points (substratum)
    Nτ::Int = 401 # [-] Number of time points
    plot_interval::Int = 400 # [-] Time points between output files
    T::Float64 # [-] End time
    L::Float64 # [-] Domain width
    H0::Float64 # [-] Initial biofilm height
    S0::Float64 = 1.0 # [-] Initial contact line position
    ε::Float64 # [-] Substratum aspect ratio
    D::Float64 # [-] Nutrient diffusion coefficient
    Q::Float64 # [-] Nutrient uptake rate
    Ψn::Float64 # [-] Biomass proliferation rate
    Ψd::Float64 # [-] Biomass death rate
    Υ::Float64 # [-] Nutrient consumption rate
    λ::Float64 # [-] Slip coefficient
end

"Data structure for experiments"
struct ExpData
    a::Float64 # [-] Agar weight percentage
    t::Vector{Float64} # [-] Time for biofilm width collection
    w::Vector{Float64} # [-] Vector of biofilm widths
    ϕ::Vector{Float64} # [-] Cell viability
    ar::Float64 # [-] Aspect ratio
end

"Apply nondimensionalisation to obtain dimensionless parameters"
function nondimensionalise(p, dp, Ds, Db)
    # Apply nondimensionalisation
    non_T = Db/dp.Xb^2*dp.T
    non_L = dp.Xs/dp.Xb
    non_H0 = dp.H0*dp.Xs/(dp.Hs*dp.Xb)
    non_ε = dp.Hs/dp.Xs
    non_D = Ds/Db
    non_Q = dp.Q*dp.Xb*dp.Xs/(dp.Hs*Db)
    return Params(T = non_T, L = non_L, H0 = non_H0, ε = non_ε, D = non_D, Q = non_Q, Ψn = p[1], Ψd = p[2], Υ = p[3], λ = p[4])
end

"Compute the objective function: distance between model and parameters for given solution"
function objective(p::Vector{T}, dp, Ds, Db, ex::ExpData) where{T}
    par = nondimensionalise(p, dp, Ds, Db)
    dist = thinfilm_slip(par, ex, false)
    return dist
end

"Function for both the objective and gradient"
function fg!(F, G, p::Vector{T}, dp, Ds, Db, ex::ExpData) where{T}
    δ::Float64 = 1e-6
    ##### Compute solution and objective function
    par = nondimensionalise(p, dp, Ds, Db)
    dist = thinfilm_slip(par, ex, false)
    ##### Compute gradient
    if G !== nothing
        # Gradient in Ψn direction
        p1 = p + δ*[1.0, 0.0, 0.0, 0.0]
        par1 = nondimensionalise(p1, dp, Ds, Db)
        dist1 = thinfilm_slip(par1, ex, false)
        G[1] = (dist1-dist)/δ
        # Gradient in Ψd direction
        p2 = p + δ*[0.0, 1.0, 0.0, 0.0] 
        par2 = nondimensionalise(p2, dp, Ds, Db)
        dist2 = thinfilm_slip(par2, ex, false)
        G[2] = (dist2-dist)/δ
        # Gradient in Υ direction
        p3 = p + δ*[0.0, 0.0, 1.0, 0.0] 
        par3 = nondimensionalise(p3, dp, Ds, Db)
        dist3 = thinfilm_slip(par3, ex, false)
        G[3] = (dist3-dist)/δ
        # Gradient in λ direction
        p4 = p + δ*[0.0, 0.0, 0.0, 1.0] 
        par4 = nondimensionalise(p4, dp, Ds, Db)
        dist4 = thinfilm_slip(par4, ex, false)
        G[4] = (dist4-dist)/δ
    end
    if F !== nothing
        return dist
    end
end

"Main function for parameter estimation"
function main()
    ##### Initialise parameters
    dp = Dimensional()
    Db = 0.24*dp.D0 # [mm^2/min] Nutrient diffusivity (biofilm)
    ##### Agar densities
    A = [0.6, 0.8, 1.2, 2.0] # Agar density
    Vp = Vector{Vector{Float64}}()
    Vf = Vector{Float64}()
    for a in A
        ##### Parameters and experimental data
        Ds = (1-0.023*a)*dp.D0 # [mm^2/min] Nutrient diffusivity (substratum)
        t, w, ϕ, ar = get_exp(a, dp, Db)
        ex = ExpData(a, t, w, ϕ, ar) # Load experimental data
        ##### Optimise parameter estimates for given experimental data
        # Global black box optimisation
        res = bboptimize(x -> objective(x, dp, Ds, Db, ex); SearchRange = [(0.0, 0.5), (0.0, 0.01), (0.0, 20.0), (0.0, 1.0)], 
            MaxFuncEvals = 500, Method = :adaptive_de_rand_1_bin)
        p0 = best_candidate(res)
        lower = [0.0, 0.0, 0.0, 0.0]
        upper = [Inf, Inf, Inf, Inf]
        # Gradient-based box-constrained optimisation
        inner_optimizer = Optim.LBFGS(m = 20, 
            linesearch = LineSearches.HagerZhang(),
            alphaguess = LineSearches.InitialHagerZhang(α0 = 1.0)) # Constructor for box-constrained solver
        @time result = Optim.optimize(Optim.only_fg!((F, G, x) -> fg!(F, G, x, dp, Ds, Db, ex)),
            lower, upper, p0, 
            Fminbox(inner_optimizer), 
            Optim.Options(outer_iterations = 3, 
                x_abstol = 1e-6, f_abstol = 1e-6, g_abstol = 1e-6, 
                outer_x_abstol = 1e-6, outer_f_abstol = 1e-6, outer_g_abstol = 1e-6, 
                show_trace = true))
        p = Optim.minimizer(result)
        f = Optim.minimum(result)
        push!(Vp, p)
        push!(Vf, f)
        writedlm("OptimalParameters.csv", Vp)
        writedlm("OptimalObjectives.csv", Vf)
        ##### Plot optimal solution
        par = nondimensionalise(p, dp, Ds, Db)
        thinfilm_slip(par, ex, true)
    end
end

@time main()