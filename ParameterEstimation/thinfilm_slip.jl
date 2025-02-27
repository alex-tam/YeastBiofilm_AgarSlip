# Solve thin-film biofilm model and compare with experiments
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

"Plot solution at a given, fixed time stamp"
function draw_solution(i, ξ, ξo, τ, S, h, ϕ, gs, gso, gb, u)
    j = τ[i]
    plot(S.*ξ, h, xlabel = L"$x$", title = "Solution (t = $j)", xlims=(0, maximum(S.*ξo)), ylims=(0, 4.0), grid = false, margin=5mm, linecolor = :black, linewidth = 2, label=L"$h$")
    plot!(S.*ξ, ϕ, linecolor = :red, linewidth = 2, label=L"$\bar{\phi}_n$")
    plot!(S.*ξ, gs, linecolor = :green, linewidth = 2, label=L"$g_s$")
    plot!(S.*ξo, gso, linecolor = :green, linewidth = 2, label=false)
    plot!(S.*ξ, gb, linecolor = :blue, linewidth = 2, label=L"$g_b$")
    plot!(S.*ξ, u, linecolor = :orange, linewidth = 2, label=L"$u$")
    savefig("sol-step-$i.pdf")
end

"Solve thin-film slip model for given model parameters"
function thinfilm_slip(par, dp, ex, output::Bool)
    ##### Summary statistics #####
    summary_statistics = Vector{Float64}()
    ### Parameters and domain ###
    ξ = range(0, 1.0, length = par.Nξ) # Computational domain
    dξ = ξ[2] - ξ[1] # Grid spacing
    τ = range(0, par.T, length = par.Nτ) # Time domain
    dτ = τ[2] - τ[1] # Time step size
    ### Initial condition ###
    S = Vector{Float64}(undef, par.Nτ)
    S[1] = 1.0 # Initial contact line position
    ξo = range(1.0, par.L/S[1], length = par.Nξ) # Initialise outer domain
    h, ϕ, gs, gso, gb = ic(par, ξ)
    u = solve_u(par, dξ, h, ϕ, gb, S[1])
    # Write initial conditions to files
    if output == true
        draw_solution(1, ξ, ξo, τ, S[1], h, ϕ, gs, gso, gb, u)
    end
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
        if S[i+1] >= S[i]
            G = LinearInterpolation(ξo_old, gso) # Create interpolations function
            ξo = range(1.0, par.L/S[i+1], length = par.Nξ) # Update outer ξ
            gso = G(ξo) # Update g_s outside biofilm
        else
            @printf("The biofilm retracted. \n"); break
        end
        # The bug in this step occurs if the biofilm retracts in the first time step, and we try to interpolate outside the domain.
        # 8. Plot solution
        if (output == true) && (mod(i, par.plot_interval) == 0)
            j = i+1
            draw_solution(j, ξ, ξo, τ, S[j], h, ϕ, gs, gso, gb, u)
        end
        idx::Int = 14/21*(par.Nτ-1)+1
        if i == idx
            idx_mid::Int = (par.Nξ+1)/2
            vol_frac = [ϕ[1], ϕ[idx_mid], ϕ[end]]
            push!(summary_statistics, norm((ex.ϕ - vol_frac)./ex.ϕ))
        end
    end
    ### Contact line position ###
    # Nondimensionalise experiment
    t_non = dp.ψn*dp.G.*ex.t # Experimental time
    w_non = ex.w./dp.Xc # Experimental half-width
    # Plot
    if output == true
        plot(τ, S, xlabel = L"$\tau$", ylabel = L"$S(\tau)$", linecolor = :black, linewidth = 2, grid = false, margin=5mm, legend = false, xlims=(0, maximum(τ)), ylims=(0, dp.Xp/dp.Xc))
        scatter!(t_non, w_non, marker =:xcross, linecolor =:black, label = "0.6% Agar")
        savefig("S.pdf")
    end
    ##### Compare with experiment
    idx_exp::Vector{Int64} = [2/21*(par.Nτ-1)+1, 4/21*(par.Nτ-1)+1, 6/21*(par.Nτ-1)+1, 8/21*(par.Nτ-1)+1, 10/21*(par.Nτ-1)+1, 12/21*(par.Nτ-1)+1, 14/21*(par.Nτ-1)+1, par.Nτ]
    S_comp = [S[i] for i in idx_exp]
    push!(summary_statistics, norm((S_comp[2:end] - w_non[2:end])./w_non[end]))
    ##### Aspect ratio
    aspect = dp.ε*maximum(h)/S[end]
    if output == true
        @printf("The aspect ratio is: %f.\n", aspect)
    end
    push!(summary_statistics, norm((aspect - ex.ar)/ex.ar))
    return sum(summary_statistics)
end