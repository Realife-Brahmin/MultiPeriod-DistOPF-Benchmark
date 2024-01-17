function plotLoadShape()
    local p2 = plot(LoadShape,
        xlabel="Hour",
        ylabel="λ",
        titlefontsize=12,
        top_margin=5mm,
        title="Load Factor λ for every hour",
        label="λ",
        xlims=(1, duration),
        xticks=1:duration,
        marker=:circle,  # Add this line for circular markers
        markercolor=:black,  # Add this line to make the markers black
        markersize=4,  # Adjust marker size as needed
        linewidth=2.5,
        minorgrid=true,
        minorgridlinestyle=:dot,
        minorgridlinewidth=2,
        minorgridcolor=:gray)

    local ext = ".png"
    local figname = "loadShape" * ext
    local fignameFull = joinpath(figureFolder, figname)
    savefig(p2, fignameFull)
end