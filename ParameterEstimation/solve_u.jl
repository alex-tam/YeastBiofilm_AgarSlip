# Solve momentum BVP for u(x,t)
# Alex Tam, 16/11/2023

function solve_u(par, dξ, h, ϕ, gb, S)
    ### Build matrix ##
    A = Matrix{Float64}(undef, par.Nξ, par.Nξ); fill!(A, 0.0) # Pre-allocate
    # Left boundary
    A[1,1] = 1.0 # Dirichlet BC (left boundary)
    # Interior grid points
    for i = 2:par.Nξ-1
        A[i,i+1] = 2/(S^2*dξ^2)*(h[i+1] + h[i])
        A[i,i] = -2/(S^2*dξ^2)*(h[i+1] + 2*h[i] + h[i-1]) - par.λ
        A[i,i-1] = 2/(S^2*dξ^2)*(h[i] + h[i-1])
    end
    # Right boundary
    A[par.Nξ, par.Nξ] = 3.0/(2*dξ) # One-sided difference
    A[par.Nξ, par.Nξ-1] = -4.0/(2*dξ) # One-sided difference
    A[par.Nξ, par.Nξ-2] = 1.0/(2*dξ) # One-sided difference
    ### Build RHS ###
    b = Vector{Float64}(undef, par.Nξ) # Pre-allocate
    # Left boundary
    b[1] = 0.0 # Dirichlet BC (left boundary)
    # Interior grid points
    for i = 2:par.Nξ-1
        b[i] = 1/(S*dξ)*(ϕ[i+1]*gb[i+1]*h[i+1] - ϕ[i-1]*gb[i-1]*h[i-1]) 
    end
    # Right boundary
    b[par.Nξ] = (S/2)*ϕ[par.Nξ]*gb[par.Nξ] # Neumann BC (right boundary)
    ### Solve linear system ###
    return A\b
end