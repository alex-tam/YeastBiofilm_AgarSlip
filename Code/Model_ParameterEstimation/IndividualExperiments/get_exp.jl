# Obtain dimensionless experimental data for thin-film model
# Alex Tam, 13/04/2025

"Obtain dimensionless experimental data for given agar density"
function get_exp(a, r, dp, Db)
    T = vec(readdlm("Results/Experimental/Width/T.csv")) # Dimensional time [hr]
    t = 60*Db/(dp.Xb^2).*T # Experimental time [-]
    ϕ = [0.8585, 0.8735, 0.9175] # Day 14 viability, converted to volume fraction
    W = vec(readdlm("Results/Experimental/Width/BiofilmHalfWidth_YPD$a-R$r.csv")) # Half-width [cm]
    w = W./(dp.Xb/10) # [-] Dimensionless half-width
    if a == 0.6 # 0.6% agar
        ar = 0.0306 # [-] Experimental aspect ratio
    elseif a == 0.8 # 0.8% agar
        ar = 0.0411 # [-] Experimental aspect ratio
    elseif a == 1.2 # 1.2% agar
        ar = 0.0577 # [-] Experimental aspect ratio
    elseif a == 2.0 # 2.0% agar
        ar = 0.0683 # [-] Experimental aspect ratio
    else
        @printf("Invalid agar density.")
    end
    return t, w, ϕ, ar
end