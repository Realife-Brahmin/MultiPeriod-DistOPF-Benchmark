# currently under implementation
using CSV
using DataFrames
include("number_to_padded_string.jl")

pv = 10
pvString = 
batt = 10
hour = 1
wd = @__DIR__

resultsFolder = joinpath(dirname(wd), "results")
configFolderName = "pv"*number_to_padded_string(pv)*"_batt"*number_to_padded_string(batt)
configFolder = joinpath(resultsFolder, configFolderName)
ext = ".csv"

filenameSummary = "summary"*string(hour)*ext
filenameVoltage = "voltages"*string(hour)*ext
filenameVoltageFull = joinpath(configFolder, filenameVoltage)
filenameSummaryFull = joinpath(configFolder, filenameSummary)
dfVoltage = CSV.read(filenameVoltageFull, DataFrame)
dfSummary = CSV.read(filenameSummaryFull, DataFrame)
dfSummary = dfSummary[end, :]
filenameBusNames = joinpath(wd, "busname.txt")
busNames = CSV.read(filenameBusNames, DataFrame, header=false)

rename!(dfVoltage, [:Bus, :BasekV, :Node1, :Magnitude1, :Angle1, :pu1]) # rename headers when I come back, remove the leading space

V_unsorted = dfVoltage.pu1
busNames