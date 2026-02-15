# Update contact line position
# Alex Tam, 16/11/2023

function solve_S(dτ, S, u, u_old)
    return S + dτ/2*(u[end] + u_old[end])
end