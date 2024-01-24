using CSV
using DataFrames
using Glob

pv = 0
pv = 10 # percentage of load buses
batt = 0
batt = 10 # percentage of load buses

function plotVoltageTimeSeries()
    # Dictionary to store DataFrames


    figureFolder = joinpath(resultsFolder, "figures")
    isdir(figureFolder) || mkdir(figureFolder)

    voltagesDict = Dict()

    configFolderName = "pv" * number_to_padded_string(pv) * "_batt" * number_to_padded_string(batt)
    configFolder = joinpath(resultsFolder, configFolderName)
    voltagesFolder = joinpath(configFolder, "voltages");

    # List all .csv files in the directory
    csv_files = glob("*.csv", voltagesFolder)

    # Iterate over the list of files
    for file in csv_files
        # Extract the unique bus (e.g., 9, 28, etc.) from the filename
        match1 = match(r"_v(\d+)_", file)
        if match1 !== nothing
            bus = match1.captures[1]

            # Read the CSV file into a DataFrame
            df = CSV.read(file, DataFrame)
            colNames = [strip(string(name)) for name in names(df)]
            rename!(df, colNames)
            # Store the DataFrame in the dictionary with a key like 'voltages9'
            voltagesDict["voltages"*bus] = df
        end
    end

    using Plots

    # Define the base voltage in kV
    kV_B = 2.4018
    kV_to_V = 1e3
    # Ensure the figure directory exists
    isdir(figureFolder) || mkdir(figureFolder)

    # Iterate over each DataFrame in the dictionary
    for (key, df) in voltagesDict
        # Extract the bus from the key (e.g., 'voltages9' -> '9')
        bus = match(r"voltages(\d+)", key).captures[1]

        # Convert voltage to per unit (pu)
        pu_values = df.V / (kV_B * kV_to_V)

        # Generate time vector (assuming evenly spaced measurements over 24 hours)
        time = range(0, stop=24, length=size(df, 1))

        # Create a plot of pu values as a function of time
        plot(time, pu_values, label="Voltage in pu", xlabel="Time (hours)", ylabel="Voltage (pu)",
            title="Voltage Time Series for bus $bus")

        # Save the plot to a file
        plot_path = joinpath(figureFolder, "voltageTimeSeries_$bus.png")
        savefig(plot_path)
    end
end

# plotVoltageTimeSeries()






