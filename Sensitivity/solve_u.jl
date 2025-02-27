# Solve momentum BVP for u(x,t)
# Alex Tam, 16/11/2023

function solve_u(par, dξ, h, ϕ, gb, S)
    ### Build matrix ##
    A = Matrix{Float64}(undef, par.Nξ, par.Nξ); fill!(A, 0.0) # Pre-allocate
    # Left boundary
    A[1,1] = 1.0 # Dirichlet BC (left boundary)
    # Interior grid points
    for i = 2:par.Nξ-1
        A[i,i+1] = 2/(S*dξ^2)*(h[i+1] + h[i])
        A[i,i] = -2/(S*dξ^2)*(h[i+1] + 2*h[i] + h[i-1]) - par.λ
        A[i,i-1] = 2/(S*dξ^2)*(h[i] + h[i-1])
    end
    # Right boundary
    A[par.Nξ, par.Nξ] = 3.0/(2*dξ) # One-sided difference
    A[par.Nξ, par.Nξ-1] = -4.0/(2*dξ) # One-sided difference
    A[par.Nξ, par.Nξ-2] = 1.0/(2*dξ) # One-sided difference
    ### Build RHS ###
    b = Vector{Float64}(undef, par.Nξ) # Pre-allocate
    Γ = surface_tension(par, dξ, h, ϕ, gb, S) # Compute surface tension
    # Left boundary
    b[1] = 0.0 # Dirichlet BC (left boundary)
    # Interior grid points
    for i = 2:par.Nξ-1
        b[i] = 2*(1+par.Ψm)/(S*dξ)*(ϕ[i+1]*gb[i+1]*h[i+1] - ϕ[i-1]*gb[i-1]*h[i-1]) - Γ[i]
    end
    # Right boundary
    b[par.Nξ] = (S/2)*(1+par.Ψm)*ϕ[par.Nξ]*gb[par.Nξ] - Γ[end] # Neumann BC (right boundary)
    ### Solve linear system ###
    return A\b
end

"Separate function for surface tension terms"
function surface_tension(par, dξ, h, ϕ, gb, S)
    Γ = Vector{Float64}(undef, par.Nξ); fill!(Γ, 0.0) # Pre-allocate
    ### Obtain surface tension ###
    # Interior grid points
    Γ[2] = par.γ*h[2]/(S^3)*(-3*h[1]+10*h[2]-12*h[3]+6*h[4]-h[5])/(2*dξ^3)
    for i = 3:par.Nξ-2
        Γ[i] = par.γ*h[i]/(S^3)*(-h[i-2]+2*h[i-1]-2*h[i+1]+h[i+2])/(2*dξ^3)
    end
    Γ[par.Nξ-1] = par.γ*h[par.Nξ-1]/(S^3)*(h[par.Nξ-4]-6*h[par.Nξ-3]+12*h[par.Nξ-2]-10*h[par.Nξ-1]+3*h[par.Nξ])/(2*dξ^3)
    # Right boundary
    Γ[end] = par.γ/(4*S)*(2*h[end]-5*h[end-1]+4*h[end-2]-h[end-3])/(dξ^2)
    return Γ
end

# Adjust the above. In my thesis, I used 6th-order accurate schemes for the surface tension, and computed higher derivatives repeating a first derivative operator.