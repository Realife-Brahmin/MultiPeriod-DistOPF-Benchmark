using CSV
using DataFrames
using Glob
using Plots

function plotVoltageTimeSeries1(resultsFolder::String)
    figureFolder = joinpath(resultsFolder, "figures")
    isdir(figureFolder) || mkdir(figureFolder)

    configs = [(0, 0), (0, 10), (10, 0), (10, 10)]
    config_colors = [:blue, :red, :green, :purple]  # Define a color for each configuration
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

                # Plot the current configuration's data on the combined plot
                plot!(combined_plot, time, pu_values, label="pv=$pv%, batt=$batt%", color=config_colors[idx])
            end
        end

        # Save the combined plot to a file
        plot_path = joinpath(figureFolder, "voltageTimeSeries_bus$bus.png")
        savefig(combined_plot, plot_path)
    end

end

# To call the function, provide the path to the results folder
plotVoltageTimeSeries1(resultsFolder)
