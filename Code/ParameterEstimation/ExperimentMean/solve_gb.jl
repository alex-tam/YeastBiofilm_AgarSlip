# Solve reaction-diffusion equations for nutrients in biofilm
# Alex Tam, 8/1/2024

"Solve for nutrient concentration in biofilm"
function solve_gb(par, dτ, dξ, ξ, h, ϕ, gs, gb, u, S)
    flag = false
    ### Build matrix ##
    A = Matrix{Float64}(undef, par.Nξ, par.Nξ); fill!(A, 0.0) # Pre-allocate
    # Left boundary
    A[1,1] = h[1] + (dτ/2)*(h[1])*(-u[3]+4*u[2]-3*u[1])/(2*S*dξ) + (dτ/2)*(h[1]+h[2])/(S^2*dξ^2) + (dτ/2)*par.Q + (dτ/2)*par.Υ*ϕ[1]*h[1]
    A[1,2] = -(dτ/2)*(h[1]+h[2])/(S^2*dξ^2)
    # Interior grid points
    for i = 2:par.Nξ-1
        A[i, i] = h[i] + (dτ/2)*( par.Q + par.Υ*ϕ[i]*h[i] + (h[i+1]+2*h[i]+h[i-1])/(2*S^2*dξ^2) + u[i]*h[i]/(S*dξ) )
        A[i, i+1] = -(dτ/2)*( h[i]*ξ[i]*u[end]/(S*dξ) + (h[i+1]+h[i])/(2*S^2*dξ^2) )
        A[i, i-1] = -(dτ/2)*( -h[i]*ξ[i]*u[end]/(S*dξ) + (h[i]+h[i-1])/(2*S^2*dξ^2) + u[i-1]*h[i-1]/(S*dξ) )
    end
    # Right boundary
    A[end, end] = 1.0
    ### Build RHS ###
    b = Vector{Float64}(undef, par.Nξ) # Pre-allocate
    # Left boundary
    b[1] = gb[1]*(h[1] - (dτ/2)*(h[1])*(-u[3]+4*u[2]-3*u[1])/(2*S*dξ) - (dτ/2)*(h[1]+h[2])/(S^2*dξ^2) - (dτ/2)*par.Q - (dτ/2)*par.Υ*ϕ[1]*h[1]) +
        gb[2]*(dτ/2)*(h[1]+h[2])/(S^2*dξ^2) + dτ*par.Q*gs[1] # Ghost node scheme
    # Interior grid points
    for i = 2:par.Nξ-1
        b[i] = h[i]*gb[i] + dτ*par.Q*gs[i] - 
            (dτ/2)*( par.Q + par.Υ*ϕ[i]*h[i] + (h[i+1]+2*h[i]+h[i-1])/(2*S^2*dξ^2) + u[i]*h[i]/(S*dξ) )*gb[i] +
            (dτ/2)*( h[i]*ξ[i]*u[end]/(2*S*dξ) + (h[i+1]+h[i])/(2*S^2*dξ^2) )*gb[i+1] +
            (dτ/2)*( -h[i]*ξ[i]*u[end]/(2*S*dξ) + (h[i]+h[i-1])/(2*S^2*dξ^2) + u[i-1]*h[i-1]/(S*dξ) )*gb[i-1]
    end
    # Right boundary
    b[end] = par.Q*gs[end]/(u[end]*(3*h[end]-4*h[end-1]+h[end-2])/(2*dξ*S) + par.Q)
    ### Solve linear system ###
    gb = A\b
    # Sanity check
    if any(gb .< 0)
        @printf("Error: Non-physical nutrient concentration in biofilm.\n")
        flag = true
        return gb, flag
    end
    return gb, flag
end