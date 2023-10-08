include("../helperFunctions.jl")

systemName = "ieee123";
numAreas = 1;
Area = 1;
folderName = "processedData/"*systemName*"/numAreas_"*string(numAreas)* "/area"*string(Area)*"/"

res = extract_values_from_files(folderName)

df = DataFrame(res)

selected_df = df[:, [:horizonTimes, :totalTimes, :nLinEqns, :nNonLinEqns, :nVars, :PLossVals]]

rename!(selected_df, :totalTimes => "totalTimes [s]", :PLossVals => "PLoss [kW]")

pretty_table(selected_df, header=names(selected_df), crop=:none)