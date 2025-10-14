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
    gso = ones(par.Nξo)
    gb = zeros(par.Nξ)
    return h, ϕ, gs, gso, gb
end

"Plot solution at a given, fixed time stamp"
function draw_solution(i, ξ, ξo, τ, S, h, ϕ, gs, gso, gb, u, a, r)
    j = τ[i]
    plot(S.*ξ, h, xlabel = L"$x$", title = "Solution (t = $j)", xlims=(0, maximum(S.*ξo)), ylims=(0, 1. + ceil(maximum(h))), grid = false, margin=5mm, linecolor = :black, linewidth = 2, label=L"$h$")
    plot!(S.*ξ, ϕ, linecolor = :red, linewidth = 2, label=L"$\phi$")
    plot!(S.*ξ, gs, linecolor = :green, linewidth = 2, label=L"$g_s$")
    plot!(S.*ξo, gso, linecolor = :green, linewidth = 2, label=false)
    plot!(S.*ξ, gb, linecolor = :blue, linewidth = 2, label=L"$g_b$")
    plot!(S.*ξ, u, linecolor = :orange, linewidth = 2, label=L"$u$")
    savefig("sol-agar-$a-R$r-step-$i.pdf")
end

"Solve thin-film slip model for given model parameters"
function thinfilm_slip(par, ex, output::Bool)
    ##### Summary statistics #####
    summary_statistics = Vector{Float64}()
    a = ex.a # Agar density for plot labelling
    r = ex.r # Replicate number for plot labelling
    idx_exp = Int.(round.(ex.t.*(par.Nτ-1)./ex.t[end])) .+ 1 # Indices for experimental comparison
    ### Parameters and domain ###
    ξ = range(0, 1.0, length = par.Nξ) # Computational domain
    dξ = ξ[2] - ξ[1] # Grid spacing
    τ = range(0, par.T, length = par.Nτ) # Time domain
    dτ = τ[2] - τ[1] # Time step size
    ### Initial condition ###
    S = Vector{Float64}(undef, par.Nτ)
    S[1] = par.S0 # Initial contact line position
    ξo = range(1.0, par.L/S[1], length = par.Nξo) # Initialise outer domain
    h, ϕ, gs, gso, gb = ic(par, ξ)
    u = solve_u(par, dξ, h, ϕ, gb, S[1]) 
    # Write initial conditions to files
    if output == true
        draw_solution(1, ξ, ξo, τ, S[1], h, ϕ, gs, gso, gb, u, a, r)
    end
    ### Time stepping ###
    for i = 1:par.Nτ-1
        u_old = u # Store old velocity for contact line
        ξo_old = ξo # Store old outer variable
        # 1. Solve height equation
        h = solve_h(par, dτ, dξ, ξ, h, ϕ, gb, u, S[i])
        if maximum(h) == h[end-1]
            @printf("Error: Spurious oscillation suspected. \n")
            return 100.0
        end
        # 2. Solve volume fraction equation
        ϕ = solve_phi(par, dτ, dξ, ξ, ϕ, gb, u, S[i])
        if any(ϕ .> 1) || any(ϕ .< 0)
            @printf("Error: Non-physical cell volume fraction. \n")
            return 100.0
        end
        # 3. Solve substratum nutrient equations (substratum)
        gs, gso, flag = solve_gs(par, dτ, dξ, ξ, ξo_old, gs, gso, gb, u, S[i])
        if flag == true
            return 100.0
        end
        # 4. Solve biofilm nutrient equation
        gb, flag = solve_gb(par, dτ, dξ, ξ, h, ϕ, gs, gb, u, S[i])
        if flag == true
            return 100.0
        end
        # 5. Solve momentum equation
        u = solve_u(par, dξ, h, ϕ, gb, S[i])
        # 6. Solve contact line equation
        S[i+1] = solve_S(dτ, S[i], u, u_old)
        # 7. Update outer domain and nutrient concentration
        x = vcat(S[i].*ξ, S[i].*ξo[2:end])
        g = vcat(gs, gso[2:end])
        if issorted(x) == false
            @printf("Error: Knot vector unsorted. \n")
            return 100.0
        end
        G = LinearInterpolation(x, g, extrapolation_bc=Interpolations.Flat()) # Create interpolations function
        ξo = range(1.0, par.L/S[i+1], length = par.Nξo) # Update outer ξ
        gs = G(S[i+1].*ξ)# Update g_s inside biofilm
        gso = G(S[i+1].*ξo) # Update g_s outside biofilm
        if S[i+1] < S[i]
            @printf("Error: The biofilm retracted. \n")
            return 100.0
        end
        if S[i+1] >= par.L
            @printf("Warning: Biofilm reached end of the domain. \n")
            return 100.0
        end
        # 8. Plot solution
        if (output == true) && (mod(i, par.plot_interval) == 0)
            j = i+1
            draw_solution(j, ξ, ξo, τ, S[j], h, ϕ, gs, gso, gb, u, a, r)
        end
        idx::Int = idx_exp[7] # Index at Day 14
        if i == idx
            idx_mid::Int = (par.Nξ+1)/2
            vol_frac = [ϕ[1], ϕ[idx_mid], ϕ[end]]
            push!(summary_statistics, norm((vol_frac - ex.ϕ)./ex.ϕ))
        end
    end
    ##### Contact line position #####
    # Plot contact line position
    if output == true
        plot(τ, S, xlabel = L"$\tau$", ylabel = L"$S(\tau)$", linecolor = :black, linewidth = 2, grid = false, margin=5mm, legend = false, xlims=(0, maximum(τ)), ylims=(0, par.L))
        scatter!(ex.t, ex.w, marker =:xcross, linecolor =:black, label = "0.6% Agar")
        savefig("S-agar-$a-R$r.pdf")
    end
    # Compare with experiment
    S_comp = [S[i] for i in idx_exp]
    push!(summary_statistics, norm((S_comp[2:end] - ex.w[2:end])./ex.w[end]))
    ##### Aspect ratio #####
    aspect = par.ε*maximum(h)/S[end]
    if output == true
        @printf("The aspect ratio is: %f.\n", aspect)
    end
    push!(summary_statistics, norm((aspect - ex.ar)/ex.ar))
    ##### Write solution to file #####
    if output == true
        writedlm("t-agar-$a-R$r.csv", τ)
        writedlm("xi-agar-$a-R$r.csv", ξ)
        writedlm("xio-agar-$a-R$r.csv", ξo)
        writedlm("h-agar-$a-R$r.csv", h)
        writedlm("phi-agar-$a-R$r.csv", ϕ)
        writedlm("gs-agar-$a-R$r.csv", gs)
        writedlm("gso-agar-$a-R$r.csv", gso)
        writedlm("gb-agar-$a-R$r.csv", gb)
        writedlm("u-agar-$a-R$r.csv", u)
        writedlm("S-agar-$a-R$r.csv", S)
    end
    ##### Compute distance #####
    return norm(summary_statistics)
end