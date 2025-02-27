# Plot solutions to thin-film biofilm model
# Alex Tam, 9/1/2024

# Load Packages
using DelimitedFiles
using Plots
using Measures
using LaTeXStrings

"Plot solution at a given, fixed time stamp"
function draw_solution(i, ξ, ξo, τ, S, h, ϕ, gs, gso, gb, u)
    j = τ[i]
    plot(S.*ξ, h, xlabel = L"$x$", title = "Solution (t = $j)", xlims=(0, maximum(S.*ξo)), ylims=(0,1.2), grid = false, margin=5mm, linecolor = :black, linewidth = 2, label=L"$h$")
    plot!(S.*ξ, ϕ, linecolor = :red, linewidth = 2, label=L"$\bar{\phi}_n$")
    plot!(S.*ξ, gs, linecolor = :green, linewidth = 2, label=L"$g_s$")
    plot!(S.*ξo, gso, linecolor = :green, linewidth = 2, label=false)
    plot!(S.*ξ, gb, linecolor = :blue, linewidth = 2, label=L"$g_b$")
    plot!(S.*ξ, u, linecolor = :orange, linewidth = 2, label=L"$u$")
    savefig("sol-step-$i.pdf")
end

"Generate all plots for one full simulation"
function draw_single(ξ, τ, plot_times)
    ### Plot configuration ###
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(fontfamily = "Computer Modern", titlefontsize = 14, guidefontsize = 18, tickfontsize = 14, legendfontsize = 12)
    ### Import solutions ###
    S = vec(readdlm("S.csv")) # Contact line position, S(τ)
    ### Plot solutions ###
    for i in plot_times
        h = vec(readdlm("h-step-$i.csv"))
        ϕ = vec(readdlm("phi-step-$i.csv"))
        gs = vec(readdlm("gs-step-$i.csv"))
        ξo = vec(readdlm("xio-step-$i.csv"))
        gso = vec(readdlm("gso-step-$i.csv"))
        gb = vec(readdlm("gb-step-$i.csv"))
        u = vec(readdlm("u-step-$i.csv"))
        draw_solution(i, ξ, ξo, τ, S[i], h, ϕ, gs, gso, gb, u)
    end
    ### Plot contact line ###
    plot(τ, S, xlabel = L"$\tau$", ylabel = L"$S(\tau)$", linecolor = :black, linewidth = 2, grid = false, margin=5mm, legend = false, xlims=(0, maximum(τ)), ylims=(0,maximum(S)))
    savefig("S.pdf")
end

"Main function to plot all required results"
function draw()
    ### Import global properties ###
    plot_times = convert(Vector{Int}, vec(readdlm("plot_times.csv"))) # Time-stamps
    ξ = vec(readdlm("xi.csv")) # Spatial variable, ξ
    τ = vec(readdlm("tau.csv")) # Time, τ
    ### Plot individual solutions ###
    draw_single(ξ, τ, plot_times)
end

@time draw()