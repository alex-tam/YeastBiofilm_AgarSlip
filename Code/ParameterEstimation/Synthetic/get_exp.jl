# Obtain dimensionless experimental data for thin-film model
# Alex Tam, 13/04/2025

"Obtain dimensionless experimental data for given agar density"
function get_exp(a, dp, Db)
    ##### Scale time
    T = vec(readdlm("Results/Experimental/Width/T.csv")) # Dimensional time [hr]
    t = 60*Db/(dp.Xb^2).*T # Experimental time [-]
    ##### Randomise
    rv = rand(Uniform(0.932, 1.05)) # Factor for random cell viability based on experimental range
    ##### Generate synthetic data
    ϕ = rv.*[0.8585, 0.8735, 0.9175] # Day 14 viability, converted to volume fraction
    w = randomise_width(a, t, dp)
    ar = randomise_aspect(a)
    return t, w, ϕ, ar
end

"Generate random aspect ratio"
function randomise_aspect(a)
    if a == 0.6
        ra = rand(Normal(1.0, 0.087)) # Factor for random aspect ratio based on cell count
        return ra*0.0306
    elseif a == 0.8
        ra = rand(Normal(1.0, 0.130)) # Factor for random aspect ratio based on cell count
        return ra*0.0411
    elseif a == 1.2
        ra = rand(Normal(1.0, 0.102)) # Factor for random aspect ratio based on cell count
        return ra*0.0577
    elseif a == 2.0
        ra = rand(Normal(1.0, 0.123)) # Factor for random aspect ratio based on cell count
        return ra*0.0683
    else
        @printf("Invalid agar density.")
    end
end

"Extract correct experimental replicates"
function get_R(y)
    if y == 0.6
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

"Generate random biofilm width"
function randomise_width(a, t, dp)
    R = get_R(a)
    W = Array{Float64}(undef, length(t), length(R))
    # Import data
    for i in eachindex(R)
        r = R[i]
        W[:, i] = vec(readdlm("Results/Experimental/Width/BiofilmHalfWidth_YPD$a-R$r.csv")) # Half-width [cm]
    end
    # Compute standard deviations
    σ = vec(std(W, dims=2))
    rw = rand(Normal())
    # Generate data
    mean = vec(readdlm("Results/Experimental/Width/BiofilmHalfWidthMean_YPD$a.csv")) # Half-width [cm]
    w = (mean .+ rw.*σ)/(dp.Xb/10) # [-] Dimensionless half-width
    return w
end