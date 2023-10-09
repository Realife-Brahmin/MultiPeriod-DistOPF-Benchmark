function plotGraphs(areaInfo, simInfo, busDataActualBusNumsTable_Area, ...
    G_Area, isRoot_Area, numAreas, CB_FullTable, numChildAreas_Area, varargin)
    
 % Default values for optional arguments
    verbose = false;
    logging = false;
    displayNetworkGraphs = true;
    displayActualBusNumbersInGraphs = false;
    systemName = "ieee123";
    loggingLocationName = "logfiles/";
    fileExtension = ".txt";
    savePlots = false;

    % Process optional arguments
    numArgs = numel(varargin);
    if mod(numArgs, 2) ~= 0
        error('Optional arguments must be specified as name-value pairs.');
    end
    
    validArgs = ["verbose", "logging", "systemName", "displayNetworkGraphs", ...
        "displayActualBusNumbersInGraphs", "savePlots"];
    
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
            case "displayNetworkGraphs"
                displayNetworkGraphs = argValue;
            case "displayActualBusNumbersInGraphs"
                displayActualBusNumbersInGraphs = argValue;
            case "savePlots"
                savePlots = argValue;
            case "fileExtension"
                fileExtension = argValue;
        end
    end
    
    macroItr = simInfo.macroItr; % completed macro-iterations, starts at 0

    saveLocationFilename = strcat(loggingLocationName , systemName, "/numAreas_", num2str(numAreas), "/optimizationLogs", fileExtension);

    fileOpenedFlag = false;

    if logging && verbose
        error("Kindly specify ONLY one of the following arguments as true: verbose and logging.")
    elseif logging && ~verbose && macroItr == 0
        fileOpenedFlag = true;
        fid = fopen(saveLocationFilename, 'w');
    elseif ~logging && macroItr == 0
        logging = verbose;
        fid = 1;
    else
        logging = false;
    end

    Area = areaInfo.Area;

    myfprintf(displayNetworkGraphs, fid,  "Plotting for Area %d:\n", Area);


    if ~isRoot_Area
        parentAreaIndices = find(CB_FullTable.childArea == Area);
        parentArea = CB_FullTable.parentArea(parentAreaIndices);
        parentNode = CB_FullTable.conBus_childAreaFrom(parentAreaIndices);
        childNode = CB_FullTable.conBus_childAreaTo(parentAreaIndices);

        myfprintf(logging, fid,  "This area's parent is Area %d.\n", parentArea);
    else

        myfprintf(logging, fid,  "This area has only the substation as its 'parent area'.\n");

    end

    if numChildAreas_Area
        childAreaIndices = find(CB_FullTable.parentArea == Area);
        childAreas = CB_FullTable.childArea(childAreaIndices);
        parentNodes = CB_FullTable.conBus_parentAreaFrom(childAreaIndices);
        childNodes = CB_FullTable.conBus_parentAreaTo(childAreaIndices);

        myfprintf(logging, fid,  "This area has %d child area", numChildAreas_Area);

        if numChildAreas_Area == 1
            myfprintf(logging, fid,  ": Area %d.\n", childAreas(numChildAreas_Area) );
        else
            myfprintf(logging, fid,  "s, which include the Areas: ");
            for i = 1:numChildAreas_Area-2
                fprintf("Area %d, ", childAreas(i) );
            end
            myfprintf(logging, fid, "Area %d and Areas %d.\n", childAreas(numChildAreas_Area-1:numChildAreas_Area) )
        end

    else
        myfprintf(logging, fid,  "This area has no child areas.\n");
    end

    if displayNetworkGraphs && macroItr == 0
        numLegendItems = 1 + ~isRoot_Area + 2*numChildAreas_Area;  
        legendList = cell(numLegendItems, 1);
        legendItr = 1;
        figureTitle = "Graph for Area " + num2str(Area) + " with " + num2str(areaInfo.N_Area) + " nodes";
        figure('Name', figureTitle);
        plotTitle = figureTitle;
        legendStr_base = strcat("Set of all branches in Area ", num2str(Area));
        legendList{legendItr} = legendStr_base;
        legendItr = legendItr + 1;
    
        if displayActualBusNumbersInGraphs
            nodeNumsForGraph = busDataActualBusNumsTable_Area.bus;
            plot1 = plot(G_Area, 'NodeLabel', nodeNumsForGraph);
        else 
            plot1 = plot(G_Area);
        end
        
        if ~isRoot_Area %POV You're the child area.
            hold on;
            scatter([NaN NaN], [NaN NaN], 'g', 'filled')
            scatter([NaN NaN], [NaN NaN], 'r', 'filled')
            hold off;
            
            highlight(plot1, parentNode, 'NodeColor', 'g');
            legendStr_parentArea_childIdx = strcat("Area ", num2str(parentArea), " to Area ", num2str(Area), ": Area ", num2str(parentArea), "'s  side");
            legendList{legendItr} = legendStr_parentArea_childIdx;
            legendItr = legendItr + 1;
            highlight(plot1, childNode, 'NodeColor', 'r');
            legendStr_childArea_childIdx = strcat("Area ", num2str(parentArea), " to Area ", num2str(Area), ": Area ", num2str(Area), "'s  side");
            legendList{legendItr} = legendStr_childArea_childIdx;
            legendItr = legendItr + 1;
        end
    
        if numChildAreas_Area %POV You're the parent area.

            hold on;
            for i = 1:numChildAreas_Area
                scatter([NaN NaN], [NaN NaN], 'm', 'filled')
                scatter([NaN NaN], [NaN NaN], 'k', 'filled')
                highlight(plot1, parentNodes(i), 'NodeColor', 'magenta');
                highlight(plot1, childNodes(i), 'NodeColor', 'black');
                legendStr_childArea_parentIdx = strcat("Area ", num2str(Area), " to Child Area(s): Area ", num2str(Area), "'s  side");
                legendList{legendItr} = legendStr_childArea_parentIdx;
                childArea = childAreas(i);
                legendStr_childArea_childIdx = strcat("Area ", num2str(Area), " to Child Area(s): Area ", num2str(childArea), "'s side");
                legendList{legendItr+1} = legendStr_childArea_childIdx;
                legendItr = legendItr + 2;
            end
            hold off;
        end

        legend(legendList, 'Location','SouthWest');
        title(plotTitle);

        if savePlots
            f = gcf;
            fileExtension = '.png';
            switch Area < 10 %single digit number
                case true
                    strArea = strcat("0", num2str(Area));
                case false
                    strArea = num2str(Area);
                otherwise
                    strArea = num2str(Area);
            end
            
            if displayActualBusNumbersInGraphs
                graphLabelType = "_fullSystemNodeNumbers";
            else
                graphLabelType = "_ownAreaNodeNumbers";
            end
            
            filenameSavedPlot = strcat("processedData/", systemName, "/numAreas_", num2str(numAreas), "/area", strArea, graphLabelType, "_network", fileExtension);
            exportgraphics(f, filenameSavedPlot, 'Resolution', 300)
        end

    end
    
    if fileOpenedFlag
        fclose(fid);
    end
    
end