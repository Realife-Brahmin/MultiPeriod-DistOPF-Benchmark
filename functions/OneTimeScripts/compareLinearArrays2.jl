using CSV
using DataFrames
using Plots
using Test

numAreas = 1
Area = rand(1:numAreas)
systemName = "ieee123"

rawDataDirB = joinpath("rawData", "$(systemName)", "numAreas_$(numAreas)", "area$(Area)")
baseDirB = joinpath("processedData", "$(systemName)", "numAreas_$(numAreas)")
println(baseDirB)
baseDir0 = baseDirB

ext = ".csv"
exttxt = ".txt"
# powerdata_string = "powerdata_no_header_csv"
powerdata_string = "powerdata"
filename_busData = joinpath(rawDataDirB, powerdata_string*exttxt)
busData = CSV.read(filename_busData, DataFrame, header=[:bus,:P_L,:Q_L,:Q_C,:P_der,:busType])
filename_branchData = joinpath(rawDataDirB, "linedata"*ext)
branchData = CSV.read(filename_branchData, DataFrame, header=true)


filePathAeq0_CSV = joinpath(baseDir0, "Aeq_0"*ext)
filePathbeq0_CSV = joinpath(baseDir0, "beq_0"*ext)
filePathlb0_CSV = joinpath(baseDir0, "lb_0"*ext)
filePathub0_CSV = joinpath(baseDir0, "ub_0"*ext)

filePathAeqB_CSV = joinpath(baseDirB, "Aeq_B"*ext)
filePathbeqB_CSV = joinpath(baseDirB, "beq_B"*ext)
filePathlbB_CSV = joinpath(baseDirB, "lb_B"*ext)
filePathubB_CSV = joinpath(baseDirB, "ub_B"*ext)

Aeq0 = CSV.read(filePathAeq0_CSV, DataFrame, header=false) |> Matrix;
AeqB = CSV.read(filePathAeqB_CSV, DataFrame, header=false) |> Matrix;

# Reading the beq_0 and beq_B vectors
beq0 = CSV.read(filePathbeq0_CSV, DataFrame, header=false)[1, :] |> Vector
beqB = CSV.read(filePathbeqB_CSV, DataFrame, header=false)[:, 1] |> Vector

# Defining the file paths for lb_0, ub_0, lb_B, and ub_B
filePathlb0_CSV = joinpath(baseDirB, "lb_0"*ext)
filePathub0_CSV = joinpath(baseDirB, "ub_0"*ext)
filePathlbB_CSV = joinpath(baseDirB, "lb_B"*ext)
filePathubB_CSV = joinpath(baseDirB, "ub_B"*ext)

# Reading the lb_0 and ub_0 vectors
lb0 = CSV.read(filePathlb0_CSV, DataFrame, header=false)[:, 1] |> Vector;
ub0 = CSV.read(filePathub0_CSV, DataFrame, header=false)[:, 1] |> Vector;

# Reading the lb_B and ub_B vectors
lbB = CSV.read(filePathlbB_CSV, DataFrame, header=false)[:, 1] |> Vector;
ubB = CSV.read(filePathubB_CSV, DataFrame, header=false)[:, 1] |> Vector;

# Plot the flipped matrix as a heatmap
hB = heatmap(AeqB, color=:blues, yflip=true, aspect_ratio=:equal, xlabel="Columns", ylabel="Rows")

# plot(hmm)
m_Area = 127
N_Area = 128
nD_Area = 85
areaInfo = Dict(:m_Area => m_Area, :N_Area => N_Area, :nD_Area => nD_Area)
@show row = rand(1:m_Area)
AeqBIndices = findall(AeqB[row, :] .!= 0)
println(index_to_variable(AeqBIndices, areaInfo))
@show AeqValues = AeqB[row, AeqBIndices]
@show beqBIdx = beqB[row]