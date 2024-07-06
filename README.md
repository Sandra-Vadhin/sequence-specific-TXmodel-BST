# sequence-specific-TXmodel-BST
This repo holds the code for generating and solving the sequence-specific transcription model of deGFP synthesis in the myTXTL cell-free system. This model was developed in Julia 1.9, but may be compatible with other versions. You will need to install Julia and the following packages:
- CSV
- DataFrames
- LinearAlgebra
- Plots
- Distributions
- DifferentialEquations
- GlobalSensitivity
- Tables
- NumericalIntegration
- Optim
- Colors
- JLD2
- BSTModelKit
- ForwardDiff
- Sundials
- DelimitedFiles
- SmoothingSplines
- Interpolations

To run the MCMC parameter estimation, type `include("TXmodel.jl") ` in the REPL.
To plot the results: `include("plotsim.jl")`.
