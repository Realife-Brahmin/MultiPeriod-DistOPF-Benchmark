using OpenDSSDirect
using CSV 
using DataFrames

# Store the current directory
original_directory = pwd()

csvfilename = "Normalized-1s-2900-pts.csv"

foldername = joinpath(".", "functions", "OneTimeScripts", "OpenDSSDirect", "5wwHs", "examples", "8500-Node")

foldername = joinpath(".", "functions", "OneTimeScripts", "OpenDSSDirect", "5wwHs", "examples", "8500-Node with some spaces in the name") # notice in the last folder that I've created a duplicate of the 8500-Node folder but inserted spaces in its name

filename = joinpath(foldername, "Master.dss")

filename_csv = joinpath(foldername, csvfilename)

# println(filename_csv)

df = CSV.File(filename_csv, header=0) |> DataFrame; # This works for both filenames, just to check that the folder with spaces in its name exists
# display(df)

dss("""
    clear
    compile $(filename)
    solve
""")

function main()
    loadnumber = Loads.First()
    kWsum = 0.0
    kvarsum = 0.0
    
    while loadnumber > 0
        kWsum += Loads.kW()
        kvarsum += Loads.kvar()
        loadnumber = Loads.Next()
    end

    kWsum, kvarsum
end

result = main()

# Reset the working directory to its original state
cd(original_directory)

println(result)
