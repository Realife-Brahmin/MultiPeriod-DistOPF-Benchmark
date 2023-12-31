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
fbus = branchData.fbus
tbus = branchData.tbus
P_L = busData.P_L
P_der = busData.P_der
Q_L = busData.Q_L
Q_C = busData.Q_C

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
kVA_B = 1000
m_Area = 127
N_Area = 128
nD_Area = 85
areaInfo = Dict(:m_Area => m_Area, :N_Area => N_Area, :nD_Area => nD_Area)

checkingRange = 1:2*m_Area+N_Area # equations which can be tested
# @show row = rand(checkingRange)
@show row = rand(2*m_Area+1:2*m_Area+N_Area)
# @show row = 382

@show iDx = getIdxFromRow(row, areaInfo)
if row <= 3*m_Area
    @show i, j = fbus[iDx], tbus[iDx]
elseif row == 2*m_Area + N_Area
    j = 1
else
    @error "floc"
end

AeqBIndices = findall(AeqB[row, :] .!= 0)
println(index_to_variable(AeqBIndices, areaInfo))
@show AeqValues = AeqB[row, AeqBIndices]
@show beqBValue = beqB[row]

if row <= m_Area
    @test beqBValue ≈ (P_L[j] - P_der[j])/kVA_B
elseif row <= 2*m_Area
    @test beqBValue ≈ (Q_L[j] - Q_C[j])/kVA_B
elseif row < 2*m_Area + N_Area
    @test beqBValue == 0
elseif row == 2*m_Area + N_Area 
    @test beqBValue == 1.03^2 
else
    @error "floc"
end