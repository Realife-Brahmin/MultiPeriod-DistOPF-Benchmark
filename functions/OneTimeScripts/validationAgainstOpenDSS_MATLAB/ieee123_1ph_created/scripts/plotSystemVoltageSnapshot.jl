function plotSystemVoltageSnapshot()
    # plot voltage profile for the t-th t
    for t = 1:duration
        local λ = LoadShape[t]
        local PSubs_kW = round(PSubs_MW_1toT[t] * MW_to_kW, digits=3)
        global PSubs[t] = PSubs_kW
        local PLosses_kW = round(PLosses_MW_1toT[t] * MW_to_kW, digits=3)
        theme(:dao)
        local nodes = 1:N
        if pv > 0
            local Irradiance = LoadShapePV[t]
            titleString="Bus Voltages for t = $(t), λ = $(λ)\n" * L"$P_{Subs}=$" * "$(PSubs_kW) kW " * L"$P_{Loss}=$" * "$(PLosses_kW) kW\n" * "with $(pv) % PVs at Irradiance $(Irradiance)" * " and $(batt) % Batteries"
        else
            titleString="Bus Voltages for t = $(t), λ = $(λ)\n" * L"$P_{Subs}=$" * "$(PSubs_kW) kW " * L"$P_{Loss}=$" * "$(PLosses_kW) kW\n" * "with $(pv) % PVs and $(batt) % Batteries"
        end

        local p1 = plot(nodes, V_1toT[nodes, t],
            xlabel="Bus Number",
            ylabel="V [pu]",
            titlefontsize=12,
            top_margin=5mm,
            title=titleString,
            label="V [pu]",
            xlims=(1, N),
            xticks=1:10:N,
            linewidth=2.5,
            marker=:circle,  # Add this line for circular markers
            markercolor=:black,  # Add this line to make the markers black
            markersize=1.0,  # Adjust marker size as needed
            minorgrid=true,
            minorgridlinestyle=:dot,
            minorgridlinewidth=2,
            minorgridcolor=:gray)

        local ext = ".png"
        local figname = "voltages$(t)" * ext
        local fignameFull = joinpath(figureFolder, figname)
        savefig(p1, fignameFull)
    end
end