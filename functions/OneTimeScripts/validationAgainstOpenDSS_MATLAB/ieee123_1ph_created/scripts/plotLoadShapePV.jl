function plotLoadShapePV()
    wd = @__DIR__
    loadShapeFolder = joinpath(dirname(wd), "data")

    ext = ".csv"
    filenameLoadShapePV = joinpath(loadShapeFolder, "LoadShape_PV" * ext)
    df_LoadShapePV = CSV.read(filenameLoadShapePV, DataFrame, header=false)
    LoadShapePV = df_LoadShapePV[:, 2]

    local p4 = plot(LoadShapePV,
        xlabel="Hour",
        ylabel=L"λ_{PV}",
        titlefontsize=12,
        top_margin=5mm,
        title="PV Irradiance for every hour",
        label=L"P_{PV}",
        xlims=(1, duration),
        xticks=1:duration,
        linewidth=2.5,
        marker=:circle,  # Add this line for circular markers
        markercolor=:black,  # Add this line to make the markers black
        markersize=4,  # Adjust marker size as needed
        minorgrid=true,
        minorgridlinestyle=:dot,
        minorgridlinewidth=2,
        minorgridcolor=:gray)

    local ext = ".png"
    local figname = "loadShapePV" * ext
    local fignameFull = joinpath(figureFolder, figname)
    savefig(p4, fignameFull)
end