# YeastBiofilm_AgarSlip
This repository contains Julia code, simulation results, and experimental data for the yeast biofilm thin-film model with slip.

## Contents

- Code
    - Model_ParameterEstimation: Julia code to obtain a set of optimal model parameters for given data.
        - ExperimentMean: Find optimal model parameters using mean experimental data for 0.6%, 0.8%, 1.2%, and 2.0% agar.
        - IndividualExperiments: Find optimal model parameters for 15 individual experimental replicates.
        - SyntheticData: Find optimal model parameters using synthetic data.
    - Model_Sensitivity: Julia code to obtain a parameter-pair heat map of the distance between the model and experiments.
    - Model_Single: Julia code to solve the thin-film mathematical model for given parameters, and plot the solution.
    - ImageProcessing: Julia code to process experimental photographs, returning .csv files and plots of the colony-biofilm half-width versus time.
    - UncertaintyPlots: Julia code the produce the box plots for the synthetic data analysis.
- Results
    - Experimental:
        - Photographs: Experimental photographs of yeast colony biofilms, labelled by agar density, replicate number, and day.
        - Width: Binary images and plots of the experimental colony-biofilm half width over time.
    - Model_Convergence: Numerical results for the optimal model parameters, for different grid spacing and time step sizes.
    - Model_ParameterEstimation: Results and plots of the optimal model parameters.
        - ExperimentalMean: Optimal parameters and plots using mean experimental data for 0.6%, 0.8%, 1.2%, and 2.0% agar.
        - IndividualExperiments: Optimal parameters and plots for 15 individual experimental replicates.
        - SyntheticData: Optimal parameters and plots for 50 synthetic experiments for each of 0.6%, 0.8%, 1.2%, and 2.0% agar.
    - Model_Sensitivity: Parameter-pair heat maps for 0.6% and 2.0% agar.
    - Model_Single: Plots and data for the optimal solutions (mean experimental data) for each of 0.6%, 0.8%, 1.2%, and 2.0% agar.

## Instructions
All code has been tested using Julia 1.12.5, which can be downloaded and installed using the instructions on the [Julia website](https://julialang.org/downloads/). To run the code in a given folder, execute the file "mail.jl".
