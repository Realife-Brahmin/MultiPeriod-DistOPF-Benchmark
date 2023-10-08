using DataFrames
using Glob
using Plots
using PrettyTables

function extract_values_from_files(folder_name)
    # Arrays to store extracted values
    totalTimes = Float64[]
    microItrs = Int[]
    horizonTimes = Int[]
    avgTimes = Float64[]
    nLinEqns = Int[]
    nNonLinEqns = Int[]
    nVars = Int[]
    objFunVals = Float64[]
    PLossVals = Float64[]
    PBattLossVals = Float64[]
    terminalSOCPenaltyVals = Float64[]

    # Locate all target files
    files = glob("Horizon_*_macroItr_*_*_optimalObjectiveFunctionValue.txt", folder_name)

    for file in files
        # Read file contents
        horizon = parse(Int, match(r"(?<=Horizon_)\d+", file).match)
        push!(horizonTimes, horizon)

        lines = readlines(file)

        # Extract values from lines
        push!(totalTimes, parse(Float64, match(r"(?<=periods took ).*?(?= \[s\])", lines[2]).match))
        push!(microItrs, parse(Int, match(r"(?<=and ).*?(?= iterations)", lines[2]).match))
        push!(avgTimes, parse(Float64, match(r"(?<=iteration: ).*?(?= \[s\])", lines[4]).match))
        push!(nLinEqns, parse(Int, match(r"(?<=Number of Linear Equations: ).*?$", lines[5]).match))
        push!(nNonLinEqns, parse(Int, match(r"(?<=Number of Nonlinear Equalities: ).*?$", lines[6]).match))
        push!(nVars, parse(Int, match(r"(?<=Number of Optimization Variables: ).*?$", lines[7]).match))
        push!(objFunVals, parse(Float64, match(r"(?<=periods = ).*?(?= \[kW\])", lines[8]).match))
        push!(PLossVals, parse(Float64, match(r"(?<=periods = ).*?(?= \[kW\])", lines[9]).match))
        push!(PBattLossVals, parse(Float64, match(r"(?<=periods = ).*?(?= \[kW\])", lines[10]).match))
        push!(terminalSOCPenaltyVals, parse(Float64, match(r"(?<=periods = ).*?(?= \[kWh\])", lines[11]).match))
        
    end

    res = (
        totalTimes = totalTimes, 
        microItrs = microItrs, 
        avgTimes = avgTimes, 
        nLinEqns = nLinEqns, 
        nNonLinEqns = nNonLinEqns, 
        nVars = nVars, 
        objFunVals = objFunVals, 
        PLossVals = PLossVals, 
        PBattLossVals = PBattLossVals, 
        terminalSOCPenaltyVals = terminalSOCPenaltyVals,
        horizonTimes = horizonTimes
    )

    return res
end

# Usage
systemName = "ieee123";
numAreas = 1;
Area = 1;
folderName = "processedData/"*systemName*"/numAreas_"*string(numAreas)* "/area"*string(Area)*"/"

res = extract_values_from_files(folderName)

df = DataFrame(res)

selected_df = df[:, [:horizonTimes, :totalTimes, :nLinEqns, :nNonLinEqns, :nVars, :PLossVals]]

# pretty_table(df, header=Symbol.(fieldnames(typeof(res))), crop=:none)
vscodedisplay(selected_df)
# pretty_table(df, crop=:none)
pretty_table(selected_df, header=names(selected_df), crop=:none)



