function plotSimulationResults(results, varargin)

    % Default values for optional arguments
        verbose = false;
        logging = false;
        systemName = "ieee123";
        fileExtension = '.png';
        saveSimulationResultPlots = true;
        processedDataFolder = strcat("processedData", filesep);
    
        % Process optional arguments
        numArgs = numel(varargin);
    
        if mod(numArgs, 2) ~= 0
            error('Optional arguments must be specified as name-value pairs.');
        end
        
        validArgs = ["verbose", "logging", "systemName", "processedDataFolder", "fileExtension", "saveSimulationResultPlots"];
        
        for i = 1:2:numArgs
            argName = varargin{i};
            argValue = varargin{i+1};
            
            if ~ischar(argName) || ~any(argName == validArgs)
                error('Invalid optional argument name.');
            end
            
            switch argName
                case "verbose"
                    verbose = argValue;
                case "logging"
                    logging = argValue;
                case "systemName"
                    systemName = argValue;
                case "processedDataFolder"
                    processedDataFolder = argValue;
                case "fileExtension"
                    fileExtension = argValue;
                case "saveSimulationResultPlots"
                    saveSimulationResultPlots = argValue;
            end
        end
         
    xSys = results.xSys;
    ySys = results.ySys;
    yArs = results.yArs;
    yNodes = results.yNodes;
    strObjectiveFunction = results.objFun;
    suffixObj = results.suffixObj;
    varyingIndependentParam = xSys.Names{1};
        
    if ~isempty(ySys)
        numTimePeriods = length(ySys.Names);
    elseif ~isempty(yArs)
        numTimePeriods = size(yArs.Vars(1), 1);
    elseif ~isempty(yNodes)
        numTimePeriods = size(yNodes.Vars(1), 1);
    else
        error("Empty results for both Systems, Areas and Nodes.");
    end
    myfprintf(logging,  "Plotter detects that %d simulations were run.\n", numTimePeriods);
    
    if ~isempty(yArs)
        numAreas = size(yArs.Vars, 2);
    elseif ~isempty(yNodes)
        numAreas = size(yArs.Vars, 2);
    else
        numAreas = 1;
    end
    myfprintf(logging,  "Plotter detects that the system had %d areas.\n", numAreas);
    
    if ~isempty(yNodes)
        NMax = size(yNodes.Vars, 3);
    else
        warning("No nodal variable?");
    end

    myfprintf(logging, "Plotter detects a maximum of %d nodes in any of the areas.\n", NMax);

    if strcmp(varyingIndependentParam, "Time Period")
        indParamString = "t";
        xLabelString = "Time-step $t$";
        titleStrAppendix = sprintf(" $t \\in [0, %d]$", numTimePeriods);
    else
        error("What is the independently varying parameter?")
    end
    
    
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
    
    
    independentVariable = xSys.Vars{1};  % Common independent variable
    titleCommon = sprintf('Objective Function: %s\n', strObjectiveFunction);  % Common part of the title string
    
    % Loop through each dependent variable and create corresponding plot
    dependentVariablesSys = ySys.Vars;
    dependentVariablesArs = yArs.Vars;

    yLabelStringsSys = ySys.yLabelNames;
    yLabelStringsArs = yArs.yLabelNames;
    legendEntriesSys = ySys.Legends;
    legendEntriesArs = yArs.yLabelNames;

    AreaNames = cell(numAreas, 1);
    for i = 1:numAreas
        AreaNames{i} = sprintf('Area %d', i);
    end

    dependentVariableNamesSys = ySys.FullNames;
    dependentVariableNamesArs = yArs.FullNames;
    figureNamesSys = ySys.FigureNames;
    figureNamesArs = yArs.FigureNames;
    
    for i = 1:numel(dependentVariablesSys)
        dependentVariable = dependentVariablesSys{i};
        figureName = figureNamesSys{i};
        dependentVariableName = dependentVariableNamesSys{i};
        legendEntrySys = dollaSign(legendEntriesSys{i});
        yLabelStrSys = dollaSign(yLabelStringsSys{i});

        figureHandle = figure('Name', strcat(dependentVariableName, " vs ", indParamString));
    
        plot(independentVariable, dependentVariable, ...
            'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
            'Color', colors{mod(i, numColours)+1}, 'LineWidth', 2.5);
    
        legend(legendEntrySys);
        titleStr = sprintf('%s%s vs %s\n', titleCommon, dependentVariableName, varyingIndependentParam) + ...
            titleStrAppendix;
        title(titleStr);
        xlabel(xLabelString);
        
        ylabel(yLabelStrSys);
        
        grid minor;
        hold off;

        saveLocation = strcat(processedDataFolder, systemName, filesep, "numAreas_", num2str(numAreas), filesep);
        filename = strcat(saveLocation, figureName, " vs ", varyingIndependentParam, '_', num2str(numTimePeriods), "_for_", suffixObj, fileExtension);
        myexportgraphics(saveSimulationResultPlots, figureHandle, filename, 'Resolution', 300);
    end

    for i = 1:numel(dependentVariablesArs)
        dependentVariable = dependentVariablesArs{i};
        figureName = figureNamesArs{i};
        dependentVariableName = dependentVariableNamesArs{i};
        legendPrefix = dollaSign(legendEntriesArs{i});
        legendEntriesAllAreas = cell(numAreas, 1);
 
        yLabelStrArs = dollaSign(yLabelStringsArs{i});

        figureHandle = figure('Name', strcat(dependentVariableName, "_vs_", indParamString));
        hold on;
        for area = 1:numAreas
            plot(independentVariable, dependentVariable(:, area), ...
                'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
                    'Color', colors{mod(area, numColours)+1}, 'LineWidth', 2.5);

            legendEntriesAllAreas{area} = strcat(legendPrefix, " ", AreaNames{area});
        end
    
        legend(legendEntriesAllAreas)
        titleStr = sprintf('%s%s vs %s\n', titleCommon, dependentVariableName, varyingIndependentParam) + ...
            titleStrAppendix;
        title(titleStr);
        xlabel(xLabelString);
        ylabel(yLabelStrArs);
        
        grid minor;
        hold off;

        saveLocation = strcat(processedDataFolder, systemName, filesep, "numAreas_", num2str(numAreas), filesep);
        filename = strcat(saveLocation, figureName, "_vs_", varyingIndependentParam, '_', num2str(numTimePeriods), "_for_", suffixObj, fileExtension);
        myexportgraphics(saveSimulationResultPlots, figureHandle, filename, 'Resolution', 300);
    end

end
