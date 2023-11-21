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

factor_W_to_kW = 1e-3

solution = Solution.Solve();
losses_kVAr = Circuit.Losses()*factor_W_to_kW
S_Substation_kVAr = Circuit.TotalPower()*factor_W_to_kW
PSubs = real(S_Substation_kVAr)
QSubs = imag(S_Substation_kVAr)
PLoss = real(losses_kVAr)
QLoss = imag(losses_kVAr)

loadnumber = Loads.First()
kWsum = 0.0
kvarsum = 0.0

while loadnumber > 0 # Loads.Next() returns 0 to signify that all loads have been taken into account
    global kWsum, kvarsum, loadnumber
    kWsum += Loads.kW()
    kvarsum += Loads.kvar()
    loadnumber = Loads.Next()
end

println("P_L = $(kWsum) kW")
println("Q_L = $(kvarsum) kVAr")
println("PLoss = $(PLoss) kW")
println("QLoss = $(QLoss) kVAr")
println("PSubs = $(PSubs) kW")
println("QSubs = $(QSubs) kVAr")