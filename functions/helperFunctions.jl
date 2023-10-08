using Glob
using Plots

function extract_values_from_files(folder_name)
    # Arrays to store extracted values
    totalTimes = Float64[]
    microItrs = Int[]
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
        lines = readlines(file)

        # Extract values from lines
        push!(totalTimes, parse(Float64, match(r"(?<=periods took ).*?(?= \[s\])", lines[2]).match))
        push!(microItrs, parse(Int, match(r"(?<=and ).*?(?= iterations)", lines[2]).match))
        # push!(avgTimes, parse(Float64, match(r"(?<=iteration: ).*?(?= \[s\])", lines[3]).match))
        # push!(nLinEqns, parse(Int, match(r"(?<=Number of Linear Equations: ).*?$", lines[4]).match))
        # push!(nNonLinEqns, parse(Int, match(r"(?<=Number of Nonlinear Equalities: ).*?$", lines[5]).match))
        # push!(nVars, parse(Int, match(r"(?<=Number of Optimization Variables: ).*?$", lines[6]).match))
        # push!(objFunVals, parse(Float64, match(r"(?<=periods = ).*?(?= \[kW\])", lines[7]).match))
        push!(PLossVals, parse(Float64, match(r"(?<=periods = ).*?(?= \[kW\])", lines[8]).match))
        push!(PBattLossVals, parse(Float64, match(r"(?<=periods = ).*?(?= \[kW\])", lines[9]).match))
        # push!(terminalSOCPenaltyVals, parse(Float64, match(r"(?<=periods = ).*?(?= \[kWh\])", lines[10]).match))
    end

    return totalTimes, microItrs, avgTimes, nLinEqns, nNonLinEqns, nVars, objFunVals, PLossVals, PBattLossVals, terminalSOCPenaltyVals
end

# Usage
systemName = "ieee123";
numAreas = 1;
Area = 1;
folderName = "processedData/"*systemName*"/numAreas_"*string(numAreas)* "/area"*string(Area)*"/"

results = extract_values_from_files(folderName)
