# include("plotSparsity.jl")
include("plot_sparsity_pattern.jl")
using CSV
using DataFrames
using Plots
using Test

# Define your variables
numAreas = 4
@show Area = rand(1:4)
systemName = "ieee123"
T = 2
# @show macroItr = rand(1:6)
@show macroItr = 1 
noBatteries = true

if !noBatteries
    battstring = "withBatteries"
else
    battstring = "withoutBatteries"
end

# Construct the relative file paths
baseDirB = joinpath("processedData", "$(systemName)", "numAreas_$(numAreas)", "area$(Area)")
baseDirA = joinpath("..","MultiPeriod-DistOPF-AlgoTesting", "processedData", "$(systemName)", "numAreas_$(numAreas)", "area$(Area)")

ext = ".csv"
filePathAeqB_CSV = joinpath(baseDirB, "Aeq_B_T_$(T)_"*battstring*"_macroItr_$(macroItr)"*ext)
filePathBeqB_CSV = joinpath(baseDirB, "beq_B_T_$(T)_"*battstring*"_macroItr_$(macroItr)"*ext)
# filePathAeqA = joinpath(baseDirA, "Aeq_A_T_$(T)_macroItr_$(macroItr)"*ext)
# filePathBeqA = joinpath(baseDirA, "beq_A_T_$(T)_macroItr_$(macroItr)"*ext)

# println(filePathAeqB)
# println(filePathAeqA)

AeqB = CSV.read(filePathAeqB_CSV, DataFrame, header=false) |> Matrix;
beqB = CSV.read(filePathBeqB_CSV, DataFrame, header=false) |> Matrix;
beqB = vec(beqB);

# AeqA = CSV.read(filePathAeqA, DataFrame, header=false) |> Matrix;
# beqA = CSV.read(filePathBeqA, DataFrame, header=false) |> Matrix;
# beqA = vec(beqA);

# plotSparsity(AeqB, beqB)

ext = ".png"
filePathAeqB = joinpath(baseDirB, "Aeq_B_T_$(T)_macroItr_$(macroItr)"*ext)

plot_sparsity_pattern(AeqB, beqB, savePlots=true, filenames=[filePathAeqB])

# filePathAeqA_in_A = joinpath(baseDirA, "Aeq_A_T_$(T)_macroItr_$(macroItr)"*ext)

# filePathAeqA_in_B = joinpath(baseDirB, "Aeq_A_T_$(T)_macroItr_$(macroItr)"*ext)

# plot_sparsity_pattern(AeqA, beqA, savePlots=true, filenames=[filePathAeqA_in_A, filePathAeqA_in_B])

# try @test AeqA ≈ AeqB atol = 1e-5
# catch
#     println("Aeq not equal.")
#     println(maximum(abs.(AeqA-AeqB)))
# end

# try @test beqA ≈ beqB atol = 1e-5
# catch
#     println("beq not equal.")
#     println(maximum(abs.(beqA-beqB)))
# end
