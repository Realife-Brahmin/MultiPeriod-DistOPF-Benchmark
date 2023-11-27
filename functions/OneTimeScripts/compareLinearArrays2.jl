using CSV
using DataFrames
using Plots
using Test

numAreas = 1
# @show Area = rand(1:4)
systemName = "ieee123"
# T = 2
# @show macroItr = rand(1:6)
# @show macroItr = 1 

baseDirB = joinpath("processedData", "$(systemName)", "numAreas_$(numAreas)")
baseDirA = baseDirB

ext = ".csv"
filePathAeq0_CSV = joinpath(baseDirB, "Aeq_0"*ext)
filePathAeqB_CSV = joinpath(baseDirB, "Aeq_B"*ext)

Aeq0 = CSV.read(filePathAeq0_CSV, DataFrame, header=false) |> Matrix; 
AeqB = CSV.read(filePathAeqB_CSV, DataFrame, header=false) |> Matrix;

mismatchLocations = Aeq0 .!= AeqB
