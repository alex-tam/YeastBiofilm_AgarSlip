# Solve mass balance equation for ϕ(x,t) using Crank-Nicolson
# Alex Tam, 8/1/2024

function solve_phi(par, dτ, dξ, ξ, ϕ, gb, u, S)
    ### Build matrix ##
    A = Matrix{Float64}(undef, par.Nξ, par.Nξ); fill!(A, 0.0) # Pre-allocate
    # Left boundary
    A[1, 1] = -3.0/(2*dξ)
    A[1, 2] = 4.0/(2*dξ)
    A[1, 3] = -1.0/(2*dξ)
    # Interior grid points
    for i = 2:par.Nξ-1
        A[i, i] = 1.0 - (dτ/2)*(par.Ψn*gb[i]*(1-ϕ[i]) - par.Ψd)
        A[i, i+1] = -dτ/(4*S*dξ)*(ξ[i]*u[end]-u[i])
        A[i, i-1] = dτ/(4*S*dξ)*(ξ[i]*u[end]-u[i])
    end
    # Right boundary
    A[end, end] = 1 - (dτ/2)*(par.Ψn*gb[end]*(1-ϕ[end]) - par.Ψd)
    ### Build RHS ###
    b = Vector{Float64}(undef, par.Nξ) # Pre-allocate
    # Left boundary
    b[1] = 0.0 # No-flux BC
    # Interior grid points
    for i = 2:par.Nξ-1
        b[i] = ϕ[i]*(1.0 + (dτ/2)*(par.Ψn*gb[i]*(1-ϕ[i]) - par.Ψd)) + dτ/(4*S*dξ)*(ξ[i]*u[end]-u[i])*(ϕ[i+1] - ϕ[i-1])
    end
    # Right boundary
    b[end] = ϕ[end]*(1.0 + (dτ/2)*(par.Ψn*gb[end]*(1-ϕ[end]) - par.Ψd))
    ### Solve linear system ###
    return A\b
end