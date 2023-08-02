function plotSimulationResults(results, varargin)

    % Default values for optional arguments
        verbose = false;
        CVR = [0; 0];
    
        systemName = "ieee123";
        fileExtension = '.png';
        saveSimulationResultPlots = true;
        processedDataFolder = strcat("processedData", filesep);
    
        % Process optional arguments
        numArgs = numel(varargin);
    
        if mod(numArgs, 2) ~= 0
            error('Optional arguments must be specified as name-value pairs.');
        end
        
        validArgs = ["verbose", "systemName", "CVR", "processedDataFolder", "fileExtension", "saveSimulationResultPlots"];
        
        for i = 1:2:numArgs
            argName = varargin{i};
            argValue = varargin{i+1};
            
            if ~ischar(argName) || ~any(argName == validArgs)
                error('Invalid optional argument name.');
            end
            
            switch argName
                case "verbose"
                    verbose = argValue;
                case "systemName"
                    systemName = argValue;
                case "CVR"
                    CVR = argValue;
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
    strObjectiveFunction = results.objFun;
    suffixObj = results.suffixObj;
    varyingIndependentParam = xSys.varyingIndependentParam;
    
    [PLoss, QLoss, delV, PSubs, QSubs, QDER, P12, Q12, V] = deal(ySys.PLoss, ySys.QLoss, ySys.delV, ySys.PSubs, ySys.QSubs, ySys.QDER, yArs.P12, yArs.Q12, yArs.V);
    
    if ~isempty(ySys)
        numSimulations = length(PLoss);
    elseif ~isempty(yArs)
        numSimulations = size(P12, 1);
    else
        error("Empty results for both Systems and Areas.");
    end
    myfprintf(verbose,  "Plotter detects that %d simulations were run.\n", numSimulations);
    
    if ~isempty(yArs)
        numAreas = size(P12, 2);
    else
        numAreas = 1;
    end
    myfprintf(verbose,  "Plotter detects that the system had %d areas.\n", numAreas);

    if strcmp(varyingIndependentParam, "Impedance")
        Impedance = xSys.Impedance;
        R = Impedance.R;
        xVals_vs_Step = R.Vals;
        mu_R = R.mu;
        sigma_R = R.sigma;
        X = Impedance.X;
        mu_X = X.mu;
        sigma_X = X.sigma;
        indParamString = "R";
        % indParamStringLatex = strcat("$", indParamString, "$");
        xLabelString = "Resistance $R \thinspace [pu]$ (Impedance Proportional)";
        titleStrAppendix = sprintf(" $R \\in [%.4f %.4f, %.4f + %.4f]$ and\n" + ...
        " $X \\in [%.5f %.5f, %.5f + %.5f]$", mu_R, -floor(numSimulations/2)*sigma_R, mu_R, floor(numSimulations/2)*sigma_R, ...
        mu_X, -floor(numSimulations/2)*sigma_X, mu_X, floor(numSimulations/2)*sigma_X);
    elseif strcmp(varyingIndependentParam, "Loading")
        lambda = xSys.lambda;
        xVals_vs_Step = lambda.Vals;
        mu_lambda = lambda.mu;
        sigma_lambda = lambda.sigma;
        mu_PL = 0.1;
        mu_QL = 0.01;
        titleStrAppendix = sprintf(" $P_L \\in [%.4f %.4f, %.4f + %.4f]$ and\n" + ...
        " $Q_L \\in [%.5f %.5f, %.5f + %.5f]$", mu_PL*mu_lambda, -floor(numSimulations/2)*sigma_lambda*mu_PL*mu_lambda, mu_PL*mu_lambda, floor(numSimulations/2)*sigma_lambda*mu_PL*mu_lambda, ...
        mu_QL*mu_lambda, -floor(numSimulations/2)*sigma_lambda*mu_QL*mu_lambda, mu_QL*mu_lambda, floor(numSimulations/2)*sigma_lambda*mu_QL*mu_lambda);
        indParamString = "lambda";
        % indParamStringLatex = strcat("$\", indParamString, "$");
        xLabelString = "Loading $\lambda$";
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
    
    
    independentVariable = xVals_vs_Step;  % Common independent variable
    titleCommon = sprintf('Objective Function: %s\n', strObjectiveFunction);  % Common part of the title string
    
    % Loop through each dependent variable and create corresponding plot
    dependentVariablesSys = {PLoss, QLoss, delV, PSubs, QSubs, QDER};
    dependentVariablesArs = {V, P12, Q12};

    yLabelStringsSys = {'$P_{Loss} \thinspace [kW]$', '$Q_{Loss} \thinspace [kVAr]$', '$\Delta V \thinspace [pu]$', "$P_{Subs} \thinspace [kW]$",  "$Q_{Subs} \thinspace [kVAr]$", "$Q_{DER} \thinspace [kVAr]$"};
    yLabelStringsArs = {'$V_1 \thinspace [pu]$', '$P_{12} \thinspace [kW]$', '$Q_{12} \thinspace [kVAr]$'};
    legendEntriesSys = {'$P_{Loss}$', '$Q_{Loss}$', '$\Delta V$', '$P_{Subs}$', "$Q_{Subs}$", '$Q_{DER}$'};
    legendEntriesArs = {"$V_1$", "$P_{12}$", "$Q_{12}$"};

    AreaNames = cell(numAreas, 1);
    for i = 1:numAreas
        AreaNames{i} = sprintf('Area %d', i);
    end

    dependentVariableNamesSys = {'Line Losses', "Line Reactive Losses", 'Voltage Deviations', 'Substation Power', "Substation Reactive Power", 'Control Variables'};
    dependentVariableNamesArs = {'Boundary Voltages', 'Boundary Real Powers', 'Boundary Reactive Powers'};
    figureNamesSys = {'lossesFigure', "reactiveLossesFigure", 'voltageDeviationsFigure', 'substationPowerFigure', "substationReactiveFigure", 'controlVariablesFigure'};
    figureNamesArs = {'boundaryVoltagesFigure', 'boundaryRealPowersFigure', 'boundaryReactivePowersFigure'};
    
    for i = 1:numel(dependentVariablesSys)
        dependentVariable = dependentVariablesSys{i};
        figureName = figureNamesSys{i};
        dependentVariableName = dependentVariableNamesSys{i};
        legendEntrySys = legendEntriesSys{i};
        yLabelStrSys = yLabelStringsSys{i};

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
        filename = strcat(saveLocation, figureName, " vs ", varyingIndependentParam, '_', num2str(numSimulations), "_for_", suffixObj, fileExtension);
        myexportgraphics(saveSimulationResultPlots, figureHandle, filename, 'Resolution', 300);
    end

    for i = 1:numel(dependentVariablesArs)
        dependentVariable = dependentVariablesArs{i};
        figureName = figureNamesArs{i};
        dependentVariableName = dependentVariableNamesArs{i};
        legendPrefix = legendEntriesArs{i};
        legendEntriesAllAreas = cell(numAreas, 1);
 
        yLabelStrArs = yLabelStringsArs{i};

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
        filename = strcat(saveLocation, figureName, "_vs_", varyingIndependentParam, '_', num2str(numSimulations), "_for_", suffixObj, fileExtension);
        myexportgraphics(saveSimulationResultPlots, figureHandle, filename, 'Resolution', 300);
    end

end
