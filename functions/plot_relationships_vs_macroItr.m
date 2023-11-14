function plot_relationships_over_time(results, simInfo, sysInfo)
    
    P12_1toT_vs_macroItr = results.P12_1toT_vs_macroItr;
    v1_1toT_vs_macroItr = results.v1_1toT_vs_macroItr;
    PLoss_allT_vs_macroItr = results.PLoss_allT_vs_macroItr;

    systemName = "ieee123";
    fileExtensionImage = '.png';
    fileExtensionData = ".csv";
    saveSimulationResultPlots = true;
    processedDataFolder = strcat("processedData", filesep);
    CBTable = sysInfo.CBTable;

    % Determine the size of the input matrix
    [numRelationships, T, totalMacroItr] = size(P12_1toT_vs_macroItr);

    % latex_interpreter
    % indParamString = "t";
    xLabelString = "Macro-iteration Number";

    numColours = 9;
    colors = cell(numColours, 1);
    
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
    alphaValues = linspace(1.0, 1.0, totalMacroItr); % From opaque to more transparent
    noBatteries = simInfo.noBatteries;
    if ~noBatteries
        battstring = "with Batteries";
    else
        battstring = "without Batteries";
    end

    numAreas = sysInfo.numAreas;
    kVA_B = sysInfo.kVA_B;
    kV_B = sysInfo.kV_B;
    % Iterate over each relationship
    for run = 1:3
        if run <= 2
            for r = 1:numRelationships
                figure; % Create a new figure for each relationship
                hold on; % Hold on to plot multiple lines
                grid minor;
        
                % Iterate over each macro iteration
                for t = 1:T
                    % Extract the data for the current relationship and macro iteration
                    if run == 2
                        data0 = P12_1toT_vs_macroItr(r, t, :);
                    elseif run == 1
                        data0 = v1_1toT_vs_macroItr(r, t, :);
                    else
                        error("Uknown thing to plot.");
                    end
    
                    data = squeeze(data0);
                    
                    if run == 2
                        dependentVariable = data*kVA_B;
                    elseif run ==1
                        % dependentVariable = data*kV_B;
                        dependentVariable = data;
                    else
                        error("Uknown thing to plot.");
                    end
                    % Plot the data over time
                    plot(1:totalMacroItr, abs(dependentVariable), ...
                    'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
                    'Color', [colors{mod(t, numColours)+1},  alphaValues(t)], 'LineWidth', 2.5, ...
                    'DisplayName', ['t = ' num2str(t)]);
    
                end
                
                parentArea = CBTable.parentArea(r);
                childArea = CBTable.childArea(r);
                % Customize the plot
                if run == 2
                    titlePre = "$P_{12}$";
                    yLabelString = "$P_{12} \, [kW]$";
                    filenamePre = "BoundaryRealPower";
                elseif run == 1
                    titlePre = "$v_{1}$";
                    yLabelString = "$v_{1} \, [pu]$";
                    filenamePre = "BoundaryVoltage";
                else
                    error("floc")
                end
    
                titleString = strcat(titlePre, " accross the horizon between Area $", num2str(CBTable.parentArea(r)),  "$ and Area $", num2str(CBTable.childArea(r)), "$ ", battstring);
                title(titleString)
                xlabel(xLabelString);
        
                % xlabel('t [units]');
                ylabel(yLabelString);
                legend show; % Show the legend
                ax = gca; % Get the current axes object
                ax.XTick = 1:totalMacroItr; % Set x-ticks to integers from 1 to totalMacroItr
        
                saveLocation = strcat(processedDataFolder, systemName, filesep, "numAreas_", num2str(numAreas), filesep);
                if ~exist(saveLocation, 'dir')
                    mkdir(saveLocation)
                end
                
                
                filenamePNG = strcat(saveLocation, filenamePre, "_vs_t_vs_macroItr_", num2str(T), "Areas_", num2str(parentArea), "_", num2str(childArea), "_", battstring,  fileExtensionImage);
                myexportgraphics(saveSimulationResultPlots, gcf, filenamePNG, 'Resolution', 300);
                filenameCSV = replace(filenamePNG, fileExtensionImage, fileExtensionData);
                writematrix(dependentVariable, filenameCSV)
        
                % saveLocation = strcat(processedDataFolder, systemName, filesep, "numAreas_", num2str(numAreas), filesep);
                % filenamePNG = strcat(saveLocation, figureNameFull, "_vs_", varyingIndependentParam, '_', num2str(T), "_for_", suffixObj, fileExtensionImage);
                % myexportgraphics(saveSimulationResultPlots, figureHandleFull, filenamePNG, 'Resolution', 300);                
            end
        end
    end


    
end
