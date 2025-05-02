# Solve mass balance equation for h(x,t) using Crank-Nicolson
# Alex Tam, 16/11/2023

function solve_h(par, dτ, dξ, ξ, h, ϕ, gb, u, S)
    ### Build matrix ##
    A = Matrix{Float64}(undef, par.Nξ, par.Nξ); fill!(A, 0.0) # Pre-allocate
    # Left boundary
    A[1,1] = -3.0/(2*dξ) # One-sided difference
    A[1,2] = 4.0/(2*dξ) # One-sided difference
    A[1,3] = -1.0/(2*dξ) # One-sided difference
    # Interior grid points
    for i = 2:par.Nξ-1
        A[i,i] = 1.0 - dτ/2*par.Ψn*ϕ[i]*gb[i] + dτ*u[i]/(2*S*dξ)
        A[i,i+1] = -dτ*u[end]*ξ[i]/(4*S*dξ) 
        A[i,i-1] = dτ*u[end]*ξ[i]/(4*S*dξ) - dτ*u[i-1]/(2*S*dξ)
    end
    # Right boundary
    A[end,end] = 1.0 - dτ/2*par.Ψn*ϕ[end]*gb[end] # One-sided Crank-Nicolson
    A[end,end-1] = dτ*u[end]/(dξ*S) - dτ*u[end-1]/(dξ*S) # One-sided Crank-Nicolson
    A[end,end-2] = -dτ*u[end]/(4*dξ*S) + dτ*u[end-2]/(4*dξ*S) # One-sided Crank-Nicolson
    ### Build RHS ###
    b = Vector{Float64}(undef, par.Nξ) # Pre-allocate
    # Left boundary
    b[1] = 0.0 # No-flux BC
    # Interior grid points
    for i = 2:par.Nξ-1
        b[i] = h[i] + dτ*u[end]*ξ[i]/(4*dξ*S)*(h[i+1]-h[i-1]) - dτ/(2*dξ*S)*(u[i]*h[i] - u[i-1]*h[i-1]) + dτ/2*par.Ψn*ϕ[i]*gb[i]*h[i]
    end
    # Right boundary
    b[end] = h[end] + dτ/2*par.Ψn*ϕ[end]*gb[end]*h[end] + dτ/(4*dξ*S)*u[end]*(3*h[end] - 4*h[end-1] + h[end-2]) - dτ/(4*dξ*S)*(3*u[end]*h[end] - 4*u[end-1]*h[end-1] + u[end-2]*h[end-2]) # One-sided Crank-Nicolson
    ### Solve linear system ###
    return A\b
end