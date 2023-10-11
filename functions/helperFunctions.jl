using DataFrames
using Glob
using Plots
using PrettyTables


function extract_values_from_files(folder_name::String, T::Int64)
    # Arrays to store extracted values
    totalTimes = Float64[]
    microItrs = Int[]
    horizonTimes = Int[]
    macroItrs = Int[]
    avgTimes = Float64[]
    nLinEqns = Int[]
    nNonLinEqns = Int[]
    nVars = Int[]
    objFunVals = Float64[]
    PLossVals = Float64[]
    PBattLossVals = Float64[]
    terminalSOCPenaltyVals = Float64[]

    # Locate all target files

    if T ≥ 1
        files = glob("Horizon_$(T)_macroItr_*_*_optimalObjectiveFunctionValue.txt", folder_name)
        # println(files)
    else
        files = glob("Horizon_\\*_macroItr_\\*_\\*_optimalObjectiveFunctionValue.txt", folder_name)
    end

    for file in files
        # Read file contents
        horizon = parse(Int, match(r"(?<=Horizon_)\d+", file).match)
        macroItr = parse(Int, match(r"(?<=macroItr_)\d+", file).match)

        push!(horizonTimes, horizon)
        push!(macroItrs, macroItr)

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
        macroItrs = macroItrs, 
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