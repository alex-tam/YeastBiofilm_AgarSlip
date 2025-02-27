# Plot solutions to thin-film biofilm model
# Alex Tam, 9/1/2024

# Load Packages
using DelimitedFiles
using Plots
using Measures
using LaTeXStrings

"Plot solution at a given, fixed time stamp"
function draw_solution(i, ξ, ξo, τ, S, h, ϕ, gs, gso, gb, u, p)
    j = τ[i]
    plot(S.*ξ, h, xlabel = L"$x$", title = "Solution (par = $p, t = $j)", xlims=(0, maximum(S.*ξo)), ylims=(0,1.2), grid = false, margin=5mm, linecolor = :black, linewidth = 2, label=L"$h$")
    plot!(S.*ξ, ϕ, linecolor = :red, linewidth = 2, label=L"$\bar{\phi}_n$")
    plot!(S.*ξ, gs, linecolor = :green, linewidth = 2, label=L"$g_s$")
    plot!(S.*ξo, gso, linecolor = :green, linewidth = 2, label=false)
    plot!(S.*ξ, gb, linecolor = :blue, linewidth = 2, label=L"$g_b$")
    plot!(S.*ξ, u, linecolor = :orange, linewidth = 2, label=L"$u$")
    savefig("sol-par-$p-step-$i.pdf")
end

"Generate all plots for one full simulation"
function draw_single(p, ξ, τ, plot_times)
    ### Plot configuration ###
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(fontfamily = "Computer Modern", titlefontsize = 14, guidefontsize = 18, tickfontsize = 14, legendfontsize = 12)
    ### Import solutions ###
    S = vec(readdlm("S-par-$p.csv")) # Contact line position, S(τ)
    ### Plot solutions ###
    for i in plot_times
        h = vec(readdlm("h-par-$p-step-$i.csv"))
        ϕ = vec(readdlm("phi-par-$p-step-$i.csv"))
        gs = vec(readdlm("gs-par-$p-step-$i.csv"))
        ξo = vec(readdlm("xio-par-$p-step-$i.csv"))
        gso = vec(readdlm("gso-par-$p-step-$i.csv"))
        gb = vec(readdlm("gb-par-$p-step-$i.csv"))
        u = vec(readdlm("u-par-$p-step-$i.csv"))
        draw_solution(i, ξ, ξo, τ, S[i], h, ϕ, gs, gso, gb, u, p)
    end
    ### Plot contact line ###
    plot(τ, S, xlabel = L"$\tau$", ylabel = L"$S(\tau)$", linecolor = :black, linewidth = 2, grid = false, margin=5mm, legend = false, xlims=(0, maximum(τ)), ylims=(0,maximum(S)))
    savefig("S-par-$p.pdf")
end

"Main function to plot all required results"
function draw()
    ### Import global properties ###
    plot_times = convert(Vector{Int}, vec(readdlm("plot_times.csv"))) # Time-stamps
    solend = plot_times[end]
    ξ = vec(readdlm("xi.csv")) # Spatial variable, ξ
    τ = vec(readdlm("tau.csv")) # Time, τ
    ### Configure parameter values ###
    parameters = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5] # Parameter values
    ContactLine = Vector{Float64}()
    Height = Vector{Vector{Float64}}()
    ### Obtain summary data and plot individual solutions ###
    for p in parameters
        draw_single(p, ξ, τ, plot_times)
        S = vec(readdlm("S-par-$p.csv")) # Contact line position, S(τ)
        h = vec(readdlm("h-par-$p-step-$solend.csv")) # Height
        push!(ContactLine, S[end])
        push!(Height, h)
    end
    ### Plot quantities versus varied parameter ###
    plot()
    for i in eachindex(parameters)
        s = ContactLine[i]
        plot!(s.*ξ, Height[i], xlabel = L"$x$", ylabel = L"$h(x, T)$", xlims=(0, 10.0), ylims=(0,1.2), grid = false, margin=5mm, linewidth = 2, label = false)
    end
    savefig("FinalProfile.pdf")
    plot(parameters, ContactLine, xlabel = L"$\lambda$", ylabel = L"$S(T)$", xlims=(0, 0.5), ylims=(0, 10.0), grid = false, margin=5mm, markershape=:xcross, markersize = 5, markercolor =:black, linecolor = :black, linewidth = 2, label = false)
    savefig("FinalSize.pdf")
end

@time draw()