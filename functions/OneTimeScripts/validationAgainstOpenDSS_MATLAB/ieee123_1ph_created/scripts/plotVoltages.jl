using CSV
using DataFrames
using LaTeXStrings
using Plots
using Plots.Measures
include("number_to_padded_string.jl")

# pv = 10 # percentage of load buses
pv = 0
# batt = 10 # percentage of load buses
batt = 0
duration = 24 # hours
N = 129 # I know this, but you can see this from the Summary or Voltage files
wd = @__DIR__

ext = ".csv"

loadShapeFolder = joinpath(dirname(wd), "data")
filenameLoadShape = joinpath(loadShapeFolder, "LoadShape1"*ext)
df_LoadShape = CSV.read(filenameLoadShape, DataFrame, header=false)
LoadShape = df_LoadShape[:, 2]
resultsFolder = joinpath(dirname(wd), "results")
configFolderName = "pv"*number_to_padded_string(pv)*"_batt"*number_to_padded_string(batt)
configFolder = joinpath(resultsFolder, configFolderName)

filenameBusNames = joinpath(wd, "busname.txt")
busNames = CSV.read(filenameBusNames, DataFrame, header=false)
busNums_OpenDSS = busNames[:, 1]
PSubs_MW_1toT, QSubs_MVAr_1toT, PLosses_MW_1toT, PLosses_pct_1toT, QLosses_MVAr_1toT = [zeros(duration) for _ ∈ 1:5]

V_1toT = zeros(N, duration)

# for hour = 1:24
for hour = 1:1
    local filenameSummary = "summary"*string(hour)*ext
    local filenameSummaryFull = joinpath(configFolder, filenameSummary)
    local dfSummary0 = CSV.read(filenameSummaryFull, DataFrame)
    local colNamesSummary = [strip(string(name)) for name in names(dfSummary0)]
    rename!(dfSummary0, colNamesSummary)
    local dfSummary = deepcopy(dfSummary0[end, :])
    # vscodedisplay(DataFrame(dfSummary))


    PSubs_MW, PSubs_MW_1toT[hour] = [dfSummary[:TotalMW] for _ ∈ 1:2]
    QSubs_MVAr, QSubs_MVAr_1toT[hour] = [dfSummary[:TotalMvar] for _ ∈ 1:2]
    PLosses_MW, PLosses_MW_1toT[hour] = [dfSummary[:MWLosses] for _ ∈ 1:2]
    PLosses_pct, PLosses_pct_1toT[hour] = [dfSummary[:pctLosses] for _ ∈ 1:2]
    QLosses_MVAr, QLosses_MVAr_1toT[hour] = [dfSummary[:MvarLosses] for _ ∈ 1:2]

    local filenameVoltage = "voltages"*string(hour)*ext
    local filenameVoltageFull = joinpath(configFolder, filenameVoltage)
    local dfVoltage = CSV.read(filenameVoltageFull, DataFrame)
    local colNamesVoltage = [strip(string(name)) for name in names(dfVoltage)]
    rename!(dfVoltage, colNamesVoltage)
    local V_unsorted = dfVoltage.pu1
    local VSubs = V_unsorted[1]
    local VGrid_unsorted = V_unsorted[2:N]
    local VGrid_sorted = VGrid_unsorted[busNums_OpenDSS]
    V_sorted = vcat(VSubs, VGrid_sorted)
    @show [V_unsorted V_sorted]
    V_1toT[:, hour] = V_sorted
end

MW_to_kW = 1000
figureFolder = joinpath(configFolder, "figures")
# plot voltage profile for the t-th hour
for t = 1:duration
    local λ = LoadShape[t]
    local PSubs_kW = PSubs_MW_1toT[t]*MW_to_kW
    local PLosses_kW = PLosses_MW_1toT[t]*MW_to_kW
    theme(:dao)
    local nodes = 1:N
    local p1 = plot(nodes, V_1toT[nodes, hour],
        xlabel="Bus Number",
        ylabel="V [pu]",
        titlefontsize=12,
        top_margin=5mm,
        title="Bus Voltages for t = $(t), λ = $(λ)\n"*L"$P_{Subs}=$"* "$(PSubs_kW) kW " * L"$P_{Loss}=$" * "$(PLosses_kW) kW\n"* "with $(pv) % PVs and $(batt) % Batteries",
        label="V [pu]",
        xlims=(1, N),
        xticks=1:10:N,
        linewidth=2.5,
        minorgrid=true,
        minorgridlinestyle=:dot,
        minorgridlinewidth=2,
        minorgridcolor=:gray)

    local ext = ".png"
    local figname = "voltages$(t)"*ext
    local fignameFull = joinpath(figureFolder, figname)
    savefig(p1, fignameFull)
end