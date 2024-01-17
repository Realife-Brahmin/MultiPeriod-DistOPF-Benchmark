function plotLoadShapeStorage()
    wd = @__DIR__
    loadShapeFolder = joinpath(dirname(wd), "data")

    ext = ".csv"
    filenameLoadShapeStorage = joinpath(loadShapeFolder, "LoadShape_Storage" * ext)
    df_LoadShapeStorage = CSV.read(filenameLoadShapeStorage, DataFrame, header=false)
    LoadShapeStorage = df_LoadShapeStorage[:, 2]

    local p3 = plot(LoadShapeStorage,
        xlabel="Hour",
        ylabel=L"P_{dischrg} - P_{chrg}",
        titlefontsize=12,
        top_margin=5mm,
        title="Battery Dispatch factor for every hour",
        label=L"P_d - P_c",
        xlims=(1, duration),
        xticks=1:duration,
        linewidth=2.5,
        marker=:circle,  # Add this line for circular markers
        markercolor=:black,  # Add this line to make the markers black
        markersize=1.0,  # Adjust marker size as needed
        minorgrid=true,
        minorgridlinestyle=:dot,
        minorgridlinewidth=2,
        minorgridcolor=:gray)

    local ext = ".png"
    local figname = "loadShapeStorage" * ext
    local fignameFull = joinpath(figureFolder, figname)
    savefig(p3, fignameFull)
end