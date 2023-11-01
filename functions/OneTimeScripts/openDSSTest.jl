using OpenDSSDirect

# filename = joinpath(dirname(dirname(OpenDSSDirect))), "examples", "8500-Node", "Master.dss") 

csvfilename = "Normalized-1s-2900-pts.csv"

# filename = joinpath("functions", "OneTimeScripts", "OpenDSSDirect", "5wwHs", "examples", "8500-Node")*"Master.dss"
filename = joinpath("functions", "OneTimeScripts", "OpenDSSDirect", "5wwHs", "examples", "8500-Node")*"Master.dss"

println(filename)

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

main()