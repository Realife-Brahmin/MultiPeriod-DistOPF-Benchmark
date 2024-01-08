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