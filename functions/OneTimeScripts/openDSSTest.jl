using OpenDSSDirect
using CSV 
using DataFrames
using Test

# Store the current directory
file_directory = @__DIR__ # directory of this .jl file

foldername = joinpath(".", "functions", "OneTimeScripts", "OpenDSSDirect", "5wwHs", "examples", "IEEE123-SinglePhase")

filename = joinpath(foldername, "prev_123JPT.dss")


dss("""
    clear
    redirect "$(filename)"
    solve
""")

solution = Solution.Solve();
losses = Circuit.Losses()
S_Substation = Circuit.TotalPower()

loadnumber = Loads.First()
kWsum = 0.0
kvarsum = 0.0

while loadnumber > 0
    global kWsum, kvarsum, loadnumber
    kWsum += Loads.kW()
    kvarsum += Loads.kvar()
    loadnumber = Loads.Next()
end

println("P_L = $kWsum")
println("Q_L = $(kvarsum)")
println("PLoss = $(real(losses))")
println("QLoss = $(imag(losses))")
println("PSubs = $(real(S_Substation))")
println("QSubs = $(imag(S_Substation))")