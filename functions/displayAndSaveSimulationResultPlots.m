function displayAndSaveSimulationResultPlots(itr, macroIterationLossesTotal, dist_timeTotal, R_max, v1_actual, saveSimulationResultPlots, Area, systemName, numAreas)
    lossesFigure = figure('Name', 'Macro-Iteration Losses');
    plot(1:itr, macroIterationLossesTotal, ...
        'Color', "#7E2F8E", 'Marker', 'o', 'MarkerFaceColor', "#7E2F8E", ...
        "MarkerEdgeColor", 'k', 'LineWidth', 2.5);
    grid minor;
    title('Macro-iteration Losses')
    xticks(1:itr)
    xticklabels(1:itr)
    xlabel('Iterations');
    ylabel('Line Loss [kW]')
    
    timesFigure = figure('Name', 'Macro-Iteration Times (varying every run)');
    plot(1:itr, dist_timeTotal, ...
        'Color', "#0072BD", 'Marker', 'o', 'MarkerFaceColor', "#0072BD", ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2.5);
    grid minor;
    title('Macro-iteration Times (varying every run)')
    xticks(1:itr)
    xticklabels(1:itr)
    xlabel('Iterations');
    ylabel('Time taken for the Iteration [s]')
    
    residualsFigure = figure('Name', 'Macro-Iteration residuals');
    plot(1:itr, R_max, ...
        'Color', "#A2142F", 'Marker', 'o', 'MarkerFaceColor', "#A2142F", ...
        "MarkerEdgeColor", 'k', 'LineWidth', 2.5);
    grid minor;
    title('Macro-iteration residuals')
    xticks(1:itr)
    xticklabels(1:itr)
    xlabel('Iterations');
    ylabel('Residuals')

    nodalVoltagesFigure = figure('Name', 'Nodal Voltages');
    plot(v1_actual, ...
        'Color', "#D95319", 'Marker', 'o', 'MarkerFaceColor', "#D95319", ...
        "MarkerEdgeColor", 'k', 'LineWidth', 2.5);
    grid minor;
    title('Nodal Voltages')
    xlabel('Nodes');
    ylabel('Voltages [pu]')

    if saveSimulationResultPlots
        fileExtension = '.png';
        switch Area < 10 %single digit number
            case true
                strArea = strcat("0", num2str(Area));
            case false
                strArea = num2str(Area);
            otherwise
                strArea = num2str(Area);
        end
        
        processedDataFolder = strcat("processedData/", systemName, "/numAreas_", num2str(numAreas), '/');
        filenameLossesFigure = strcat(processedDataFolder, 'macroItrLosses', fileExtension);
        filenameTimesFigure = strcat(processedDataFolder, 'macroItrTimes', fileExtension);
        filenameResidualsFigure = strcat(processedDataFolder, 'R_max', fileExtension);
        filenameNodalVoltagesFigure = strcat(processedDataFolder, 'v1_actual', fileExtension);
        exportgraphics(lossesFigure, filenameLossesFigure, 'Resolution', 300)
        exportgraphics(timesFigure, filenameTimesFigure, 'Resolution', 300)
        exportgraphics(residualsFigure, filenameResidualsFigure, 'Resolution', 300)
        exportgraphics(nodalVoltagesFigure, filenameNodalVoltagesFigure, 'Resolution', 300)
    end
end