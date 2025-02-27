# Solve thin-film biofilm model with slip using front-fixing method
# Alex Tam, 07/11/2023

"Initial conditions"
function ic(par, ξ)
    h = Vector{}(undef, par.Nξ)
    for i in eachindex(ξ)
        h[i] = par.H0*(1-ξ[i]^2)
    end
    ϕ = ones(par.Nξ)
    gs = ones(par.Nξ)
    gso = ones(par.Nξ)
    gb = zeros(par.Nξ)
    return h, ϕ, gs, gso, gb
end

"Solve thin-film slip model for given model parameters"
function thinfilm_slip(par)
    ### Parameters and domain ###
    ξ = range(0, 1.0, length = par.Nξ); dξ = ξ[2] - ξ[1] # Computational domain
    τ = range(0, par.T, length = par.Nτ); dτ = τ[2] - τ[1] # Time domain
    writedlm("xi.csv", ξ); writedlm("tau.csv", τ) # Write domain to files
    plot_times = Vector{Int}() # Vector of time-steps at which to store results
    push!(plot_times, 1)
    ### Initial condition ###
    S = Vector{Float64}(undef, par.Nτ); S[1] = 1.0 # Initial contact line position
    ξo = range(1.0, par.L/S[1], length = par.Nξ) # Initialise outer domain
    h, ϕ, gs, gso, gb = ic(par, ξ)
    u = solve_u(par, dξ, h, ϕ, gb, S[1])
    # Write initial conditions to files
    writedlm("h-step-1.csv", h); 
    writedlm("phi-step-1.csv", ϕ); 
    writedlm("gs-step-1.csv", gs); 
    writedlm("xio-step-1.csv", ξo); 
    writedlm("gso-step-1.csv", gso); 
    writedlm("gb-step-1.csv", gb); 
    writedlm("u-step-1.csv", u)
    ### Time stepping ###
    for i = 1:par.Nτ-1
        u_old = u # Store old velocity for contact line
        ξo_old = ξo # Store old outer variable
        # 1. Solve height equation
        h = solve_h(par, dτ, dξ, ξ, h, ϕ, gb, u, S[i])
        # 2. Solve volume fraction equation
        ϕ = solve_phi(par, dτ, dξ, ξ, ϕ, gb, u, S[i])
        # 3. Solve substratum nutrient equations (substratum)
        gs, gso = solve_gs(par, dτ, dξ, ξ, ξo_old, gs, gso, gb, u, S[i])
        # 4. Solve biofilm nutrient equation
        gb = solve_gb(par, dτ, dξ, ξ, h, ϕ, gs, gb, u, S[i])
        # 5. Solve momentum equation
        u = solve_u(par, dξ, h, ϕ, gb, S[i])
        # 6. Solve contact line equation
        S[i+1] = solve_S(dτ, S[i], u, u_old)
        # 7. Update outer domain and nutrient concentration
        G = LinearInterpolation(ξo_old, gso) # Create interpolations function
        ξo = range(1.0, par.L/S[i+1], length = par.Nξ) # Update outer ξ
        gso = G(ξo) # Update g_s outside biofilm
        # 8. Write solution to files
        if mod(i, par.plot_interval) == 0
            j = i+1
            writedlm("h-step-$j.csv", h)
            writedlm("phi-step-$j.csv", ϕ)
            writedlm("gs-step-$j.csv", gs)
            writedlm("xio-step-$j.csv", ξo)
            writedlm("gso-step-$j.csv", gso)
            writedlm("gb-step-$j.csv", gb)
            writedlm("u-step-$j.csv", u)
            push!(plot_times, j)
            writedlm("plot_times.csv", plot_times)
        end
    end
    writedlm("S.csv", S)
end