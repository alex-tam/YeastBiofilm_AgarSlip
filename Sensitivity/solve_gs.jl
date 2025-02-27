# Solve reaction-diffusion equations for nutrients in substratum
# Alex Tam, 8/1/2024

"Solve for nutrient concentration in substratum (Newton's method)"
function solve_gs(par, dτ, dξ, ξ, ξo, gs, gso, gb, u, S)
    dξo = ξo[2] - ξo[1] # Outer grid spacing
    # Newton iteration to match concentration at ξ = 1
    gr = gs[end] # Initial guess of nutrient concentration at ξ = 1
    for iterations = 1:10
        gr_old = gr # Store initial guess
        gs_test = solve_gs_inner(par, dτ, dξ, gr, ξ, gs, gb, u, S)
        gso_test = solve_gs_outer(par, dτ, dξo, gr, ξo, gso, u, S)
        gs_pert = solve_gs_inner(par, dτ, dξ, gr + par.ε, ξ, gs, gb, u, S)
        gso_pert = solve_gs_outer(par, dτ, dξo, gr + par.ε, ξo, gso, u, S)
        # Apply continuity condition
        f = (-3*gso_test[1] + 4*gso_test[2] - gso_test[3])/(2*dξo) -  (3*gs_test[end] - 4*gs_test[end-1] + gs_test[end-2])/(2*dξ)
        fp = (-3*gso_pert[1] + 4*gso_pert[2] - gso_pert[3])/(2*dξo) -  (3*gs_pert[end] - 4*gs_pert[end-1] + gs_pert[end-2])/(2*dξ)
        dfdr = (fp - f)/par.ε
        # Update guess
        gr = gr - f/dfdr
        # Check tolerance
        if gr - gr_old < par.ε
            break
        end
        if iterations == 10
            @printf("Nutrient matching did not converge. \n")
        end
    end
    # Sanity check
    if gr > 1
        @printf("Nutrient matching error occurred.\n")
    end
    # Obtain solutions with corrected gr
    gs = solve_gs_inner(par, dτ, dξ, gr, ξ, gs, gb, u, S)
    gso = solve_gs_outer(par, dτ, dξo, gr, ξo, gso, u, S)
    return gs, gso
end

"Solve for nutrient concentration on ξ ∈ [0, 1]"
function solve_gs_inner(par, dτ, dξ, gr, ξ, gs, gb, u, S)
    ### Build matrix ##
    A = Matrix{Float64}(undef, par.Nξ, par.Nξ); fill!(A, 0.0) # Pre-allocate
    # Left boundary
    A[1, 1] = 1.0 + (dτ/2)*( 2*par.D/(S^2*dξ^2) + par.D*par.Qs )
    A[1, 2] = -(dτ/2)*( 2*par.D/(S^2*dξ^2) )
    # Interior grid points
    for i = 2:par.Nξ-1
        A[i, i] = 1.0 + (dτ/2)*( 2*par.D/(S^2*dξ^2) + par.D*par.Qs )
        A[i, i+1] = -(dτ/2)*( ξ[i]*u[end]/(2*S*dξ) + par.D/(S^2*dξ^2) )
        A[i, i-1] = -(dτ/2)*( -ξ[i]*u[end]/(2*S*dξ) + par.D/(S^2*dξ^2) )
    end
    # Right boundary
    A[end, end] = 1.0 # Dirichlet condition
    ### Build RHS ###
    b = Vector{Float64}(undef, par.Nξ) # Pre-allocate
    # Left boundary
    b[1] = gs[1] + dτ*par.D*par.Qs*gb[1] - (dτ/2)*( 2*par.D/(S^2*dξ^2) + par.D*par.Qs )*gs[1] + (dτ/2)*( 2*par.D/(S^2*dξ^2) )*gs[2] # No-flux BC
    # Interior grid points
    for i = 2:par.Nξ-1
        b[i] = gs[i] + dτ*par.D*par.Qs*gb[i] - (dτ/2)*( 2*par.D/(S^2*dξ^2) + par.D*par.Qs )*gs[i] + (dτ/2)*( ξ[i]*u[end]/(2*S*dξ) + par.D/(S^2*dξ^2) )*gs[i+1] + (dτ/2)*( -ξ[i]*u[end]/(2*S*dξ) + par.D/(S^2*dξ^2) )*gs[i-1]
    end
    # Right boundary
    b[end] = gr # Dirichlet condition
    ### Solve linear system ###
    return A\b
end

"Solve for nutrient concentration on ξ ∈ [1, ξ_max]"
function solve_gs_outer(par, dτ, dξ, gr, ξ, gs, u, S)
    ### Build matrix ##
    A = Matrix{Float64}(undef, par.Nξ, par.Nξ); fill!(A, 0.0) # Pre-allocate
    # Left boundary
    A[1, 1] = 1.0
    # Interior grid points
    for i = 2:par.Nξ-1
        A[i, i] = 1.0 + (dτ/2)*( 2*par.D/(S^2*dξ^2) )
        A[i, i+1] = -(dτ/2)*( ξ[i]*u[end]/(2*S*dξ) + par.D/(S^2*dξ^2) )
        A[i, i-1] = -(dτ/2)*( -ξ[i]*u[end]/(2*S*dξ) + par.D/(S^2*dξ^2) )
    end
    # Right boundary
    A[end, end] = 1.0 + (dτ/2)*( 2*par.D/(S^2*dξ^2) )
    A[end, end-1] = -(dτ/2)*( 2*par.D/(S^2*dξ^2) )
    ### Build RHS ###
    b = Vector{Float64}(undef, par.Nξ) # Pre-allocate
    # Left boundary
    b[1] = gr
    # Interior grid points
    for i = 2:par.Nξ-1
        b[i] = gs[i] - (dτ/2)*( 2*par.D/(S^2*dξ^2) )*gs[i] + (dτ/2)*( ξ[i]*u[end]/(2*S*dξ) + par.D/(S^2*dξ^2) )*gs[i+1] + (dτ/2)*( -ξ[i]*u[end]/(2*S*dξ) + par.D/(S^2*dξ^2) )*gs[i-1]
    end
    # Right boundary
    b[end] = gs[end] - (dτ/2)*( 2*par.D/(S^2*dξ^2) )*gs[end] + (dτ/2)*( 2*par.D/(S^2*dξ^2) )*gs[end-1]
    ### Solve linear system ###
    return A\b
end