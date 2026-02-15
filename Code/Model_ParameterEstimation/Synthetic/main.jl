# Parameter estimation for thin-film slip model
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
using Distributions

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
    D0::Float64 = 4.04e-2 # [mm^2/min] Glucose diffusivity in water
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
    Ψn::Float64 # [-] Biomass proliferation rate
    Ψd::Float64 # [-] Biomass death rate
    D::Float64 # [-] Nutrient diffusion coefficient
    Q::Float64 # [-] Nutrient uptake rate
    Υ::Float64 # [-] Nutrient consumption rate
    λ::Float64 # [-] Slip coefficient
end

"Data structure for experimental data"
struct ExpData
    a::Float64 # [-] Agar weight percentage
    r::Int # [-] Replicate number
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
    return Params(T = non_T, L = non_L, H0 = non_H0, ε = non_ε, Ψn = p[1], Ψd = p[2], D = non_D, Q = p[3], Υ = p[4], λ = p[5])
end

"Compute the objective function: distance between model and experiment for given parameters"
function objective(p::Vector{T}, dp, Ds, Db, ex::ExpData) where{T}
    par = nondimensionalise(p, dp, Ds, Db)
    dist = thinfilm_slip(par, ex, false)
    return dist
end

"Function to compute both the objective function and gradient"
function fg!(F, G, p::Vector{T}, dp, Ds, Db, ex::ExpData) where{T}
    δ::Float64 = 1e-6
    ##### Compute solution and objective function
    par = nondimensionalise(p, dp, Ds, Db)
    dist = thinfilm_slip(par, ex, false)
    ##### Compute gradient
    if G !== nothing
        # Ψn component
        p1 = p + δ*[1.0, 0.0, 0.0, 0.0, 0.0]
        par1 = nondimensionalise(p1, dp, Ds, Db)
        dist1 = thinfilm_slip(par1, ex, false)
        G[1] = (dist1-dist)/δ
        # Ψd component
        p2 = p + δ*[0.0, 1.0, 0.0, 0.0, 0.0] 
        par2 = nondimensionalise(p2, dp, Ds, Db)
        dist2 = thinfilm_slip(par2, ex, false)
        G[2] = (dist2-dist)/δ
        # Q component
        p3 = p + δ*[0.0, 0.0, 1.0, 0.0, 0.0] 
        par3 = nondimensionalise(p3, dp, Ds, Db)
        dist3 = thinfilm_slip(par3, ex, false)
        G[3] = (dist3-dist)/δ
        # Υ component
        p4 = p + δ*[0.0, 0.0, 0.0, 1.0, 0.0] 
        par4 = nondimensionalise(p4, dp, Ds, Db)
        dist4 = thinfilm_slip(par4, ex, false)
        G[4] = (dist4-dist)/δ
        # λ component
        p5 = p + δ*[0.0, 0.0, 0.0, 0.0, 1.0] 
        par5 = nondimensionalise(p5, dp, Ds, Db)
        dist5 = thinfilm_slip(par5, ex, false)
        G[5] = (dist5-dist)/δ
    end
    if F !== nothing
        return dist
    end
end

"Main function for parameter estimation"
function main()
    ##### Establish plots and parameters
    gr() # Load GR plotting back-end
    default(titlefont = (14, "Computer Modern"), guidefont = (18, "Computer Modern"), tickfont = (14, "Computer Modern")) # Plot settings
    dp = Dimensional()
    Db = 0.24*dp.D0 # [mm^2/min] Nutrient diffusivity (biofilm)
    ##### Agar densities
    A = [0.6, 0.8, 1.2, 2.0] # Agar density
    for a in A
        Vp = Vector{Vector{Float64}}()
        Vf = Vector{Float64}()
        W = Vector{Vector{Float64}}()
        Phi = Vector{Vector{Float64}}()
        AR = Vector{Float64}()
        for r = 1:25
            ##### Parameters and synthetic data
            Ds = (1-0.023*a)*dp.D0 # [mm^2/min] Nutrient diffusivity (substratum)
            t, w, ϕ, ar = get_exp(a, dp, Db)
            ex = ExpData(a, r, t, w, ϕ, ar) # Generate synthetic data
            push!(W, w*(dp.Xb/10)) # Store synthetic data
            push!(Phi, ϕ) # Store synthetic data
            push!(AR, ar) # Store synthetic data
            writedlm("Data_W-$a.csv", W)
            writedlm("Data_Phi-$a.csv", Phi)
            writedlm("Data_AR-$a.csv", AR)
            ##### Optimise parameter estimates for given data
            # Global black-box optimisation
            res = bboptimize(x -> objective(x, dp, Ds, Db, ex); SearchRange = [(0.0, 1.0), (0.0, 0.02), (0.0, 10.0), (0.0, 10.0), (0.0, 5.0)], 
                MaxFuncEvals = 100, Method = :adaptive_de_rand_1_bin)
            p0 = best_candidate(res)
            lower = [0.0, 0.0, 0.0, 0.0, 0.0]
            upper = [Inf, Inf, Inf, Inf, Inf]
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
            writedlm("OptimalParameters-$a.csv", Vp)
            writedlm("OptimalObjectives-$a.csv", Vf)
        end
    end
end

@time main()