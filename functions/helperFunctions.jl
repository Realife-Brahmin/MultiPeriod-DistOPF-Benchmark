using DataFrames
using Glob
using Plots
using PrettyTables

"""
    extract_values_from_files(folder_name::String)

Extract various simulation results from text files contained within a specified folder.

# Arguments
- `folder_name::String`: The path to the folder containing the target `.txt` files.

# Returns
- `NamedTuple`: A named tuple containing the following fields:
  * `totalTimes`: An array of total times for the optimization in each file.
  * `microItrs`: An array of the number of iterations taken in each file.
  * `avgTimes`: An array of the average times per iteration for each file.
  * `nLinEqns`: An array of the number of linear equations in each file.
  * `nNonLinEqns`: An array of the number of nonlinear equalities in each file.
  * `nVars`: An array of the number of optimization variables in each file.
  * `objFunVals`: An array of the objective function values for each file.
  * `PLossVals`: An array of the real power line loss values for each file.
  * `PBattLossVals`: An array of the battery power loss values for each file.
  * `terminalSOCPenaltyVals`: An array of the average SOC level constraint violation values for each file.
  * `horizonTimes`: An array indicating the time horizon for each file.

# Notes
The function searches for files matching the pattern `Horizon\_*_macroItr\_*_*\_optimalObjectiveFunctionValue.txt`
in the specified folder and extracts the above values from them.

# Example
```julia
results = extract_values_from_files("path/to/folder")
"""
function extract_values_from_files(folder_name::String)
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