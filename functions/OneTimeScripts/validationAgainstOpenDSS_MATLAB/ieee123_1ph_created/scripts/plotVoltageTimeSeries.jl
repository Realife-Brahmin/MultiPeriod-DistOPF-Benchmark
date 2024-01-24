using Colors
using CSV
using DataFrames
using Glob
using Plots

function plotVoltageTimeSeries(resultsFolder::String)
    figureFolder = joinpath(resultsFolder, "figures")
    isdir(figureFolder) || mkdir(figureFolder)

    configs = [(0, 0), (0, 10), (10, 0), (10, 10)]
    # config_colors = [:blue, :red, :green, :purple]  # Define a color for each configuration
    # config_colors = [:gray, :blue, :orange, :green]
    config_colors = [:darkgray, :royalblue, :orangered, :forestgreen]
    # config_linestyles = [:solid, :dot, :dash, :dashdot]
    config_linestyles = [:solid, :solid, :solid, :solid]
    config_markers = [(:hexagon, 6, 0.6), (:circle, 6, 0.6), (:star, 9, 0.6), (:square, 6, 0.6)]
    config_markerColors = [:black, :black, :black, :black]
    config_linewidths = [5.0, 4.0, 3.0, 4.0]
    kV_B = 2.4018
    kV_to_V = 1e3

    # Collect all unique bus numbers from the first configuration to determine which plots to make
    firstConfigFolderName = "pv$(number_to_padded_string(configs[1][1]))_batt$(number_to_padded_string(configs[1][2]))"
    firstConfigFolder = joinpath(resultsFolder, firstConfigFolderName)
    firstVoltagesFolder = joinpath(firstConfigFolder, "voltages")
    bus_files = glob("*.csv", firstVoltagesFolder)
    bus_numbers = unique(match(r"_v(\d+)_", file).captures[1] for file in bus_files if match(r"_v(\d+)_", file) !== nothing)

    # ... (previous code remains unchanged)

    # Now iterate over each bus number to generate the plots
    for bus in bus_numbers
        # Initialize a plot with a title and labels, but no data
        combined_plot = plot(title="Voltage Time Series for Bus $bus", xlabel="Time (hours)", ylabel="Voltage (pu)", legend=:topright)

        for (idx, (pv, batt)) in enumerate(configs)
            configFolderName = "pv$(number_to_padded_string(pv))_batt$(number_to_padded_string(batt))"
            configFolder = joinpath(resultsFolder, configFolderName)
            voltagesFolder = joinpath(configFolder, "voltages")
            csv_files = glob("*_v$(bus)_*.csv", voltagesFolder)

            # Check if there is a file for the current bus in this configuration
            if !isempty(csv_files)
                df = CSV.read(csv_files[1], DataFrame)
                colNames = [strip(string(name)) for name in names(df)]
                rename!(df, colNames)

                # Convert voltage to per unit (pu)
                pu_values = df.V / (kV_B * kV_to_V)

                # Generate time vector (assuming evenly spaced measurements over 24 hours)
                time = range(0, stop=24, length=size(df, 1))

                Plots.theme(:dao)
                # Plot the current configuration's data on the combined plot
                yl = "V_{$(bus)}"
                plot!(combined_plot, 
                time, 
                pu_values, 
                xlabel="Hour",
                ylabel = L"V_{%$(bus)} \, [pu]",
                titlefontsize=12,
                top_margin=4mm,
                title="Voltage Time-series for Bus $(bus)\n"*"for Various Grid Edge Device Configurations",
                label= L"V \; [pu]" * " for pv=$pv%, batt=$batt%",
                legend=:bottomleft,
                xlims=(1, duration),
                xticks=1:duration,
                linewidth=config_linewidths[idx],
                linestyle=config_linestyles[idx],
                color=config_colors[idx],
                marker=config_markers[idx],  # Add this line for circular markers
                markevery=5,
                markercolor=config_markerColors[idx],  # Add this line to make the markers black
                # markersize=4,  # Adjust marker size as needed
                minorgrid=true,
                minorgridlinestyle=:dot,
                minorgridlinewidth=1,
                minorgridcolor=:gray
                )
            end
        end

        # Save the combined plot to a file
        plot_path = joinpath(figureFolder, "voltageTimeSeries_bus$bus.png")
        savefig(combined_plot, plot_path)
    end

end

# # To call the function, provide the path to the results folder
# plotVoltageTimeSeries(resultsFolder)
