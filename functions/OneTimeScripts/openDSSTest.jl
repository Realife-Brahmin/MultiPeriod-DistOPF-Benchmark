using OpenDSSDirect
using CSV 
using DataFrames
using Test

# Store the current directory
file_directory = @__DIR__ # directory of this .jl file
wd = dirname(dirname(file_directory)) # goes to the root directory of my workspace
println(wd)

foldername = joinpath(".", "functions", "OneTimeScripts", "OpenDSSDirect", "5wwHs", "examples", "IEEE123-SinglePhase")

filename = joinpath(foldername, "prev_123JPT.dss")

filename_csv = joinpath(foldername, csvfilename)


dss("""
    clear
    redirect "$(filename)"
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


# println(result)
