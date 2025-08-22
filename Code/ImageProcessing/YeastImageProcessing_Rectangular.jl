# Process images and compute biofilm width for rectangular yeast experiments
# Alex Tam, 15/04/2024

# Packages
using Printf # C-style printing macro
using FileIO # File handling tools
using Images # Image processing tools
using Plots # Plotting library
using LaTeXStrings # LaTeX in plot labels
using Measures # Allow adjustment of plot margins
using StatsBase # Statistical tools
using DelimitedFiles

"Generate file names for images"
function filename(y, r, d)
    fw = string("Projects_Current/Biofilm_ThinFilmSlip/Experiments/YPDMat_Width/YPD", y, "_R", r, "_D", d, ".jpeg")
    fs = string("Projects_Current/Biofilm_ThinFilmSlip/Experiments/YPDMat_Scale/YPD", y, "_R", r, "_D", d, ".jpeg")
    return fw, fs
end

"Extract correct experimental replicates"
function get_R(y)
    if y == 0.4
        return [1, 2, 3, 4]
    elseif y == 0.6
        return [1, 2, 3, 4]
    elseif y == 0.8
        return [1, 2, 4, 5]
    elseif y == 1.2
        return [1, 2, 3, 4]
    elseif y == 2.0
        return [1, 2, 3]
    else
        @printf("Invalid agar density.")
    end
end

"Compute biofilm width (in pixels)"
function compute_width(im)
    L = label_components(im) # Label connected regions in the binary image
    count = maximum(component_lengths(L)[1:end]) # Count total number of pixels in the biofilm (component length vector indexes from zero)
    height = size(im)[1] # Get height in pixels from size (height, width)
    return count/height # Mean width (pixels)
end

"Compute scaling factor from biofilm image"
function compute_scale(im)
    return size(im)[1]
end

"Obtain averaged data from multiple replicates"
function yeast_av(S_all, T)
    m = Vector{Float64}() # Mean biofilm width for all replicates of the same agar density
    ϵ⁻ = Vector{Float64}() # Lower bound for all replicates of the same agar density
    ϵ⁺ = Vector{Float64}() # Upper bound for all replicates of the same agar density
    for d in eachindex(T)
        v = Vector{Float64}()
        for i in eachindex(S_all)
            push!(v, S_all[i][d])
        end
        push!(m, mean(v))
        push!(ϵ⁻, mean(v) - minimum(v))
        push!(ϵ⁺, maximum(v) - mean(v))
    end
    return m, ϵ⁻, ϵ⁺
end

"Main function for image processing operations"
function main()
    ##### Plot settings #####
    gr(); plot() # Load GR plotting back-end
    default(titlefont = (12, "Computer Modern"), guidefont = (14, "Computer Modern"), tickfont = (12, "Computer Modern"), legendfont = (12, "Computer Modern")) # Plot settings
    ##### Manually enter experiment properties #####
    YPD = [0.6, 0.8, 1.2, 2.0] # Agar density
    D = [2, 4, 6, 8, 10, 12, 14, 21] # Time of photograph (day)
    T = [43., 95., 138., 186., 234.5, 289., 339.5, 513.] # Time of photograph (hour)
    writedlm("T.csv", T)
    ##### Pre-allocate #####
    M = Vector{}() # Mean widths for all agar densities and time
    E⁻ = Vector{}() # Lower bounds for all agar densities and time
    E⁺ = Vector{}() # Upper bounds for all agar densities and time
    ##### Process individual images #####
    for y in YPD # Loop over agar density
        S_all = Vector{Vector{Float64}}() # Initialise empty vector of vectors for multiple replicates
        R = get_R(y)
        for r in R # Loop over replicates
            S = Vector{Float64}() # Empty vector of biofilm half-width width for a single replicate (vs. time)
            for d in D # Loop over photographs of same replicate
                file_width, file_scale = filename(y, r, d) # Obtain file names for photographs
                im = Gray.(load(file_width)) # Load image (width) in grayscale using FileIO
                imw = binarize(im, Otsu()) # Binarize image using Otsu's method
                save("Binary_YPD$y-R$r-D$d.pdf", imw) # Save binary image
                ims = load(file_scale) # Load image (scale)
                w = compute_width(imw) # Compute biofilm width (pixels)
                p = compute_scale(ims) # Compute length scale (pixels per 10cm)
                push!(S, 0.5*10*w/p) # Store half-width
            end
            writedlm("BiofilmHalfWidth_YPD$y-R$r.csv", S)
            push!(S_all, S)
            # Generate plot for single experiment
            scatter(T, S, 
                marker=:xcross, markerstrokewidth=2, markersize = 5, color=:black, linewidth=2,
                xlabel = L"t \; [\textrm{hr}]", ylabel = L"S(t) \; [\textrm{cm}]", label=false, 
                xlims = (0, 520), ylims=(0, 2.0), margins = 5mm)
            savefig("BiofilmHalfWidth_YPD$y-R$r.pdf")
        end
        ##### Generate data for all replicates of the same agar density #####
        m, ϵ⁻, ϵ⁺ = yeast_av(S_all, T)
        writedlm("BiofilmHalfWidthMean_YPD$y.csv", m)
        ##### Generate plots for all replicates of the same agar density #####
        scatter(T, m, yerr=(ϵ⁻, ϵ⁺), 
            marker=:xcross, markerstrokewidth=2, markersize = 5, color=:black, linewidth=2,
            xlabel = L"t \; [\textrm{hr}]", ylabel = L"S(t) \; [\textrm{cm}]", label=false, 
            xlims = (0, 520), ylims=(0, 2.0), margins = 5mm)
        savefig("BiofilmHalfWidth_MeanScatter_YPD$y.pdf")
        plot(T, m, ribbon=(ϵ⁻, ϵ⁺),
            fillalpha = 0.5, marker=:xcross, markerstrokewidth=2, markersize = 5, color=:black, linewidth=2,
            xlabel = L"t \; [\textrm{hr}]", ylabel = L"S(t) \; [\textrm{cm}]", label=false, 
            xlims = (0, 520), ylims=(0, 2.0), margins = 5mm)
        savefig("BiofilmHalfWidth_MeanRibbon_YPD$y.pdf")
        ##### Store global data #####
        push!(M, m)
        push!(E⁻, ϵ⁻)
        push!(E⁺, ϵ⁺)
    end
    # Plot results for all time, replicates, and agar densities on the same axes
    for i in eachindex(M)
        m = M[i]
        ϵ⁻ = E⁻[i]
        ϵ⁺ = E⁺[i]
        y = YPD[i]
        if i == 1
            plot(T, m, yerr=(ϵ⁻, ϵ⁺), 
                marker=:xcross, markerstrokewidth=2, markersize = 5, markerstrokecolor=:auto, lc=:auto, linewidth=2,
                xlabel = L"t \; [\textrm{hr}]", ylabel = L"S(t) \; [\textrm{cm}]", label= "$y% Agar", 
                xlims = (0, 520), ylims=(0, 2.0), margins = 5mm)
        else
            plot!(T, m, yerr=(ϵ⁻, ϵ⁺), marker=:xcross, markerstrokewidth=1, markersize = 3, markerstrokecolor=:auto, lc=:auto, linewidth=2, label= "$y% Agar")
        end
    end
    savefig("BiofilmHalfWidth_MeanScatter_All.pdf")
end

@time main()