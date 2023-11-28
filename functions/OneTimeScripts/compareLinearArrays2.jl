using CSV
using DataFrames
using Plots
using Test

numAreas = 1
systemName = "ieee123"


baseDirB = joinpath("processedData", "$(systemName)", "numAreas_$(numAreas)")
baseDir0 = baseDirB

ext = ".csv"
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


misAeq = Aeq0 .!= AeqB

@show nMMAeq = sum(misAeq) # number of MisMatches
mismatch_indices = findall(misAeq)
# mismatch_indices_ij = [(i[1], i[2]) for i in mismatch_indices]
mismatch_indices_ij = [(i[2], i[1]) for i in mismatch_indices]

# Flip the matrix vertically
misAeq_flipped = reverse(misAeq, dims=1)

# Plot the flipped matrix as a heatmap
hmm = heatmap(misAeq_flipped, color=:blues, yflip=true, aspect_ratio=:equal, xlabel="Columns", ylabel="Rows")

h0 = heatmap(misAeq_flipped, color=:blues, yflip=true, aspect_ratio=:equal, xlabel="Columns", ylabel="Rows")

plot(hmm)

misbeq = beq0 .!= beqB
@show sum(misbeq)
mislb  = lb0 .!= lbB
@show sum(mislb)
misub = ub0 .!= ubB
@show sum(misub)