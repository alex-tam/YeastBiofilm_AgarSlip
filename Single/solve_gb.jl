# Solve reaction-diffusion equations for nutrients in biofilm
# Alex Tam, 8/1/2024

"Solve for nutrient concentration in biofilm"
function solve_gb(par, dτ, dξ, ξ, h, ϕ, gs, gb, u, S)
    ### Build matrix ##
    A = Matrix{Float64}(undef, par.Nξ, par.Nξ); fill!(A, 0.0) # Pre-allocate
    # Left boundary
    A[1,1] = -3.0/(2*dξ) 
    A[1,2] = 4.0/(2*dξ)
    A[1,3] = -1.0/(2*dξ)
    # Interior grid points
    for i = 2:par.Nξ-1
        A[i, i] = par.Pe*h[i] + (dτ/2)*( par.Qb + par.Υ*ϕ[i]*h[i] + (h[i+1]+2*h[i]+h[i-1])/(2*S^2*dξ^2) + par.Pe*u[i]*h[i]/(S*dξ) )
        A[i, i+1] = -(dτ/2)*( par.Pe*h[i]*ξ[i]*u[end]/(2*S*dξ) + (h[i+1]+h[i])/(2*S^2*dξ^2) )
        A[i, i-1] = -(dτ/2)*( -par.Pe*h[i]*ξ[i]*u[end]/(2*S*dξ) + (h[i]+h[i-1])/(2*S^2*dξ^2) + par.Pe*u[i-1]*h[i-1]/(S*dξ) )
    end
    # Right boundary
    A[end, end] = 3.0/(2*dξ) # One-sided difference
    A[end, end-1] = -4.0/(2*dξ) # One-sided difference
    A[end, end-2] = 1.0/(2*dξ) # One-sided difference
    ### Build RHS ###
    b = Vector{Float64}(undef, par.Nξ) # Pre-allocate
    # Left boundary
    b[1] = 0.0 # No-flux
    # Interior grid points
    for i = 2:par.Nξ-1
        b[i] = par.Pe*h[i]*gb[i] + dτ*par.Qb*gs[i] - 
            (dτ/2)*( par.Qb + par.Υ*ϕ[i]*h[i] + (h[i+1]+2*h[i]+h[i-1])/(2*S^2*dξ^2) + par.Pe*u[i]*h[i]/(S*dξ) )*gb[i] +
            (dτ/2)*( par.Pe*h[i]*ξ[i]*u[end]/(2*S*dξ) + (h[i+1]+h[i])/(2*S^2*dξ^2) )*gb[i+1] +
            (dτ/2)*( -par.Pe*h[i]*ξ[i]*u[end]/(2*S*dξ) + (h[i]+h[i-1])/(2*S^2*dξ^2) + par.Pe*u[i-1]*h[i-1]/(S*dξ) )*gb[i-1]
    end
    # Right boundary
    b[end] = 0.0 # No-flux
    ### Solve linear system ###
    return A\b
end