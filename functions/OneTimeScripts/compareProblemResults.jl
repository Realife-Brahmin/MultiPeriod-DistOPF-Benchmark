include("../helperFunctions.jl")

systemName = "ieee123";
T = 1;

numAreas = 4;

for Area = 1:numAreas
    # Area = 1;
    if Area > numAreas
        @warn "Area number $(Area) ≥ numArea $(numAreas)"
    end
    local folderName = "processedData/"*systemName*"/numAreas_"*string(numAreas)* "/area"*string(Area)*"/"
    println(folderName)
    local res = extract_values_from_files(folderName, T)

    local df = DataFrame(res)

    local selected_df = df[:, [:horizonTimes, :macroItrs, :totalTimes, :nLinEqns, :nNonLinEqns, :nVars, :PLossVals]]

    rename!(selected_df, :totalTimes => "totalTimes [s]", :PLossVals => "PLoss [kW]")

    sort!(selected_df, order(:macroItrs))

    pretty_table(selected_df, header=names(selected_df), crop=:none)
end