using Pkg
Pkg.activate("./openDSSFilesCreation")
using CSV
using DataFrames
using Glob
using LaTeXStrings
using Plots
using Plots.Measures
include("number_to_padded_string.jl")
include("plotLoadShape.jl")
include("plotLoadShapePV.jl")
include("plotLoadShapeStorage.jl")
include("plotPSubs.jl")
include("plotSystemVoltageSnapshot.jl")
include("plotVoltageTimeSeries.jl")

pv = 0
# pv = 10 # percentage of load buses
batt = 0
batt = 10 # percentage of load buses

duration = 24 # hours
N = 129 # I know this, but you can see this from the Summary or Voltage files
wd = @__DIR__

ext = ".csv"

loadShapeFolder = joinpath(dirname(wd), "data")
filenameLoadShape = joinpath(loadShapeFolder, "LoadShape1" * ext)
df_LoadShape = CSV.read(filenameLoadShape, DataFrame, header=false)
LoadShape = df_LoadShape[:, 2]

resultsFolder = joinpath(dirname(wd), "results")
configFolderName = "pv" * number_to_padded_string(pv) * "_batt" * number_to_padded_string(batt)
configFolder = joinpath(resultsFolder, configFolderName)
figureFolder = joinpath(configFolder, "figures")

voltagesFolder = joinpath(configFolder, "voltages");


plotLoadShape()

if batt > 0
    plotLoadShapeStorage()
end

if pv > 0
    filenameLoadShapePV = joinpath(loadShapeFolder, "LoadShape_PV" * ext)
    df_LoadShapePV = CSV.read(filenameLoadShapePV, DataFrame, header=false)
    LoadShapePV = df_LoadShapePV[:, 2]
    plotLoadShapePV()
end

filenameBusNames = joinpath(wd, "busname.txt")
busNames = CSV.read(filenameBusNames, DataFrame, header=false)
busNums_OpenDSS = busNames[:, 1]
PSubs_MW_1toT, QSubs_MVAr_1toT, PLosses_MW_1toT, PLosses_pct_1toT, QLosses_MVAr_1toT = [zeros(duration) for _ ∈ 1:5]

V_1toT, V_unsorted_1toT = [zeros(N, duration) for _ ∈ 1:2]

for t = 1:duration
    local filenameSummary = "summary" * string(t) * ext
    local filenameSummaryFull = joinpath(configFolder, filenameSummary)
    local dfSummary0 = CSV.read(filenameSummaryFull, DataFrame)
    local colNamesSummary = [strip(string(name)) for name in names(dfSummary0)]
    rename!(dfSummary0, colNamesSummary)
    local dfSummary = deepcopy(dfSummary0[end, :])
    # vscodedisplay(DataFrame(dfSummary))


    PSubs_MW, PSubs_MW_1toT[t] = [dfSummary[:TotalMW] for _ ∈ 1:2]
    QSubs_MVAr, QSubs_MVAr_1toT[t] = [dfSummary[:TotalMvar] for _ ∈ 1:2]
    PLosses_MW, PLosses_MW_1toT[t] = [dfSummary[:MWLosses] for _ ∈ 1:2]
    PLosses_pct, PLosses_pct_1toT[t] = [dfSummary[:pctLosses] for _ ∈ 1:2]
    QLosses_MVAr, QLosses_MVAr_1toT[t] = [dfSummary[:MvarLosses] for _ ∈ 1:2]

    local filenameVoltage = "voltages" * string(t) * ext
    local filenameVoltageFull = joinpath(configFolder, filenameVoltage)
    local dfVoltage = CSV.read(filenameVoltageFull, DataFrame)
    local colNamesVoltage = [strip(string(name)) for name in names(dfVoltage)]
    rename!(dfVoltage, colNamesVoltage)
    local V_unsorted = dfVoltage.pu1
    V_unsorted_1toT[:, t] = V_unsorted
    local VSubs = V_unsorted[1]
    local VGrid_unsorted = V_unsorted[2:N]
    # Create an array of pairs (bus number, voltage)
    local bus_voltage_pairs = [(busNums_OpenDSS[i], VGrid_unsorted[i]) for i in 1:N-1]
    # Sort the pairs by bus number
    local sorted_pairs = sort(bus_voltage_pairs, by=x->x[1])

    # Extract the sorted voltages
    local VGrid_sorted = [pair[2] for pair in sorted_pairs]


    # local VGrid_sorted = VGrid_unsorted[busNums_OpenDSS]
    local V_sorted = vcat(VSubs, VGrid_sorted)
    # @show [V_unsorted V_sorted]
    global V_1toT[:, t] = V_sorted
end

PSubs = zeros(duration)
MW_to_kW = 1000

plotSystemVoltageSnapshot()

plotPSubs()

# figureFolderCombined = joinpath(resultsFolder, "figures")
plotVoltageTimeSeries(resultsFolder)