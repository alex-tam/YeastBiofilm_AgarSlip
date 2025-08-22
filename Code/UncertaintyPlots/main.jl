# Plot convergence and uncertainty results
# Alex Tam, 22/07/2025

# Import packages
using DelimitedFiles
using StatsPlots
using LaTeXStrings
using Measures

"Main function for generating results"
function main()
    ##### Establish plots
    gr() # Load GR plotting back-end
    default(titlefont = (14, "Computer Modern"), guidefont = (18, "Computer Modern"), tickfont = (14, "Computer Modern")) # Plot settings
    # Mean optimal parameters
    v1 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/ExperimentalMean/v1/OptimalParameters.csv")
    v2 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/ExperimentalMean/v2/OptimalParameters.csv")
    v3 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/ExperimentalMean/v3/OptimalParameters.csv")
    v4 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/ExperimentalMean/v4/OptimalParameters.csv")
    v5 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/ExperimentalMean/v5/OptimalParameters.csv")
    mean = (v1 .+ v2 .+ v3 .+ v4 .+ v5)./5
    writedlm("MeanOptimalParameters.csv", mean)
    v1 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/ExperimentalMean/v1/OptimalObjectives.csv")
    v2 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/ExperimentalMean/v2/OptimalObjectives.csv")
    v3 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/ExperimentalMean/v3/OptimalObjectives.csv")
    v4 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/ExperimentalMean/v4/OptimalObjectives.csv")
    v5 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/ExperimentalMean/v5/OptimalObjectives.csv")
    mean = (v1 .+ v2 .+ v3 .+ v4 .+ v5)./5
    writedlm("MeanOptimalObjectives.csv", mean)
    # Synthetic Data Box Plots
    a1 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/SyntheticData/OptimalParameters-0.6.csv")
    a2 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/SyntheticData/OptimalParameters-0.8.csv")
    a3 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/SyntheticData/OptimalParameters-1.2.csv")
    a4 = readdlm("Projects_Current/Biofilm_ThinFilmSlip/Results/ParameterEstimation/SyntheticData/OptimalParameters-2.0.csv")
    boxplot(["0.6"], a1[:,1], fillalpha=0.75, linewidth=2, ylims = (0.0, 1.0), legend=false, xlabel = "Agar Density (%)", ylabel=L"$\Psi_n$", bar_width = 0.5, margins=3mm, size=(300, 500))
    boxplot!(["0.8"], a2[:,1], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    boxplot!(["1.2"], a3[:,1], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    boxplot!(["2.0"], a4[:,1], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    savefig("SyntheticBox_Psin.pdf")
    boxplot(["0.6"], a1[:,2], fillalpha=0.75, linewidth=2, ylims = (0.0, 0.015), legend=false, xlabel = "Agar Density (%)", ylabel=L"$\Psi_d$", bar_width = 0.5, margins=3mm, size=(300, 500))
    boxplot!(["0.8"], a2[:,2], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    boxplot!(["1.2"], a3[:,2], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    boxplot!(["2.0"], a4[:,2], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    savefig("SyntheticBox_Psid.pdf")
    boxplot(["0.6"], a1[:,3], fillalpha=0.75, linewidth=2, ylims = (0.0, 25.0), legend=false, xlabel = "Agar Density (%)", ylabel=L"$Q^{*}$", bar_width = 0.5, margins=3mm, size=(300, 500))
    boxplot!(["0.8"], a2[:,3], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    boxplot!(["1.2"], a3[:,3], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    boxplot!(["2.0"], a4[:,3], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    savefig("SyntheticBox_Q.pdf")
    boxplot(["0.6"], a1[:,4], fillalpha=0.75, linewidth=2, ylims = (0.0, 12.0), legend=false, xlabel = "Agar Density (%)", ylabel=L"$\Upsilon$", bar_width = 0.5, margins=3mm, size=(300, 500))
    boxplot!(["0.8"], a2[:,4], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    boxplot!(["1.2"], a3[:,4], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    boxplot!(["2.0"], a4[:,4], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    savefig("SyntheticBox_Upsilon.pdf")
    boxplot(["0.6"], a1[:,5], fillalpha=0.75, linewidth=2, ylims = (0.0, 6.0), legend=false, xlabel = "Agar Density (%)", ylabel=L"$\lambda^{*}$", bar_width = 0.5, margins=3mm, size=(300, 500))
    boxplot!(["0.8"], a2[:,5], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    boxplot!(["1.2"], a3[:,5], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    boxplot!(["2.0"], a4[:,5], fillalpha=0.75, linewidth=2, bar_width = 0.5)
    savefig("SyntheticBox_Lambda.pdf")
end

@time main()