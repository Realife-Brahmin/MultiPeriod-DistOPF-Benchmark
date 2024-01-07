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

filenameVoltage = "voltages"*string(hour)*ext
filename = joinpath(configFolder, filenameVoltage)

dfVoltage = CSV.read(filename, DataFrame)

filenameBusNames = joinpath(wd, "busname.txt")
busNames = CSV.read(filenameBusNames, DataFrame, header=false)
