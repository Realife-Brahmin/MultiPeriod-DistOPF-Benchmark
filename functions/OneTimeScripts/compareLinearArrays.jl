include("plotSparsity.jl")
include("plot_sparsity_pattern.jl")
using CSV
using DataFrames

# Define your variables
numAreas = 4
Area = 2
systemName = "ieee123"
T = 1
macroItr = 1

# Construct the relative file paths
baseDirB = joinpath("processedData", "$(systemName)", "numAreas_$(numAreas)", "area$(Area)")
baseDirA = joinpath("..","MultiPeriod-DistOPF-AlgoTesting", "processedData", "$(systemName)", "numAreas_$(numAreas)", "area$(Area)")

filePathAeqB = joinpath(baseDirB, "Aeq_B_T_$(T)_macroItr_$(macroItr).csv")
filePathBeqB = joinpath(baseDirB, "beq_B_T_$(T)_macroItr_$(macroItr).csv")
filePathAeqA = joinpath(baseDirA, "Aeq_A_T_$(T)_macroItr_$(macroItr).csv")
filePathBeqA = joinpath(baseDirA, "beq_A_T_$(T)_macroItr_$(macroItr).csv")

# println(filePathAeqB)
# println(filePathAeqA)

AeqB = CSV.read(filePathAeqB, DataFrame, header=false) |> Matrix;
beqB = CSV.read(filePathBeqB, DataFrame, header=false) |> Matrix;
beqB = vec(beqB);

AeqA = CSV.read(filePathAeqA, DataFrame, header=false) |> Matrix;
beqA = CSV.read(filePathBeqA, DataFrame, header=false) |> Matrix;
beqA = vec(beqA);

# plotSparsity(AeqB, beqB)

plot_sparsity_pattern(AeqB, beqB)

