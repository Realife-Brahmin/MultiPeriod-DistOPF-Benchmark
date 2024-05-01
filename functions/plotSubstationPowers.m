function plotSubstationPowers(results)
    
    simInfo = results.simInfo;
    sysInfo = results.sysInfo;

    systemName = "ieee123";
    simNatureString = simInfo.simNatureString;
    fileExtensionImage = '.png';
    fileExtensionData = ".csv";
    saveSimulationResultPlots = true;
    processedDataFolder = strcat("processedData", filesep);

    latex_interpreter
    
    suffixObj = simInfo.suffixObj;
    
    numColours = 9;
    colors = cell(numColours, 1);
    T = simInfo.T;
    totalMacroItr = simInfo.macroItr+1;

    colors{1} = [0.929, 0.694, 0.125];   % Yellow
    colors{2} = [0.494, 0.184, 0.556];   % Purple
    colors{3} = [0.466, 0.674, 0.188];   % Green
    colors{4} = [0.301, 0.745, 0.933];   % Light Blue
    colors{5} = [0.635, 0.078, 0.184];   % Dark Red
    colors{6} = [0.850, 0.325, 0.098];   % Orange
    colors{7} = [0.929, 0.694, 0.498];   % Tan
    colors{8} = [0.494, 0.694, 0.796];   % Blue
    colors{9} = [0.929, 0.612, 0.329];   % Salmon
    % Define some colors for the lines
    % colors = lines(totalMacroItr); % MATLAB's built-in colormap for lines
    alphaValues = linspace(1.0, 1.0, max(totalMacroItr, T)); % From opaque to more transparent
    % noBatteries = simInfo.noBatteries;
    DER_percent = simInfo.DER_percent;
    Batt_percent = simInfo.Batt_percent;

    battstring = simInfo.battstring;
    battstringTitle = strcat("with $", num2str(DER_percent), "\%$ PVs and $", num2str(Batt_percent), "\%$ Batteries");


    numAreas = sysInfo.numAreas;
    kVA_B = sysInfo.kVA_B;
    kV_B = sysInfo.kV_B;
    % Iterate over each relationship

    xLabelString = "Time Period $t$";

    figure; % Create a new figure for each relationship
    hold on; % Hold on to plot multiple lines
    grid minor;
    
    PSubs_1toT = results.PSubs_1toT;
    data0 = PSubs_1toT(1:T);

    data = squeeze(data0);

    dependentVariable = data*kVA_B;

    % Plot the data over time
    plot(1:T, abs(dependentVariable), ...
        'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
        'Color', [colors{mod(totalMacroItr, numColours)+1},  alphaValues(totalMacroItr)], 'LineWidth', 2.5, ...
        'DisplayName', 'Converged values');

    % Customize the plot

    titlePre = "$P_{Subs}$";
    yLabelString = "$P_{Subs} \, [kW]$";
    filenamePre = "SubstationRealPowers";


    titleString = [
        strcat(titlePre, " across the horizon for the system using"), 
        strcat(simNatureString, " ", battstringTitle) 
    ];
    title(titleString)
    xlabel(xLabelString);

    % xlabel('t [units]');
    ylabel(yLabelString);
    legend show; % Show the legend
    ax = gca; % Get the current axes object
    ax.XTick = 1:T; % Set x-ticks to integers from 1 to totalMacroItr

    saveLocation = strcat(processedDataFolder, systemName, filesep, "numAreas_", num2str(numAreas), filesep);
    if ~exist(saveLocation, 'dir')
        mkdir(saveLocation)
    end


    filenamePNG = strcat(saveLocation, filenamePre, "_vs_t_", num2str(T), "_", suffixObj, "_", battstring,  fileExtensionImage);
    myexportgraphics(saveSimulationResultPlots, gcf, filenamePNG, 'Resolution', 300);
    filenameCSV = replace(filenamePNG, fileExtensionImage, fileExtensionData);
    writematrix(dependentVariable, filenameCSV)

end