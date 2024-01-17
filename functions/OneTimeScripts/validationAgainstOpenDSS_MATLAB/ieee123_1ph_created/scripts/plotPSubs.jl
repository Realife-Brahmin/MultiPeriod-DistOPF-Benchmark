function plotPSubs()
    local p5 = plot(PSubs,
        xlabel="Hour",
        ylabel="P_{Subs}",
        titlefontsize=12,
        top_margin=5mm,
        title="Power Borrowed from the Substation for every hour",
        label="P_{Subs}",
        xlims=(1, duration),
        xticks=1:duration,
        linewidth=2.5,
        minorgrid=true,
        minorgridlinestyle=:dot,
        minorgridlinewidth=2,
        minorgridcolor=:gray)

    local ext = ".png"
    local figname = "PSubs" * ext
    local fignameFull = joinpath(figureFolder, figname)
    savefig(p5, fignameFull)
end