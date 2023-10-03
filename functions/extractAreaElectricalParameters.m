function [busDataTable_pu_Area, branchDataTable_Area, ...
    edgeMatrix_Area, R_Area, X_Area] = ...
    extractAreaElectricalParameters(Area, t, macroItr, isRoot_Area, systemName, numAreas, CB_FullTable, numChildAreas_Area, varargin)

 % Default values for optional arguments
    verbose = false;
    logging = false;
    displayTables = false;
    displayNetworkGraphs = true;
    displayActualBusNumbersInGraphs = false;
    saveLocationName = "logfiles/";
    fileExtension = ".txt";
    savePlots = true;
    kV_B = 4.16/sqrt(3);
    kVA_B = 1000;
    load_mult = 1;
    gen_mult = 1;
    systemDataFolder = strcat("rawData", filesep, systemName, filesep, "numAreas_", num2str(numAreas), filesep);

    % Process optional arguments
    numArgs = numel(varargin);
    if mod(numArgs, 2) ~= 0
        error('Optional arguments must be specified as name-value pairs.');
    end
    
    validArgs = ["verbose", "logging", "displayTables", "displayNetworkGraphs", ...
        "displayActualNumbersInGraphs", "savePlots", "kV_B", ...
        "kVA_B", "load_mult", "gen_mult"];
    
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
            case "displayTables"
                displayTables = argValue;
            case "displayNetworkGraphs"
                displayNetworkGraphs = argValue;
            case "displayActualNumbersInGraphs"
                displayActualBusNumbersInGraphs = argValue;
            case "savePlots"
                savePlots = argValue;
            case "kV_B"
                kV_B = argValue;
            case "kVA_B"
                kVA_B = argValue;
            case "load_mult"
                load_mult = argValue;
            case "gen_mult"
                gen_mult = argValue;
            case "systemDataFolder"
                systemDataFolder = argValue;
            case "fileExtension"
                fileExtension = argValue;
        end
    end
    
    saveLocationFilename = strcat(saveLocationName , systemName, "/numAreas_", num2str(numAreas), "/optimizationLogs", fileExtension);

    fileOpenedFlag = false;

    if logging && verbose
        error("Kindly specify ONLY one of the following arguments as true: verbose and logging.")
    elseif logging && ~verbose && macroItr == 1
        fileOpenedFlag = true;
        if t == 1
            fid = fopen(saveLocationFilename, 'w');
        else
            fid = fopen(saveLocationFilename, 'a');
        end
    elseif ~logging && macroItr == 1
        logging = verbose;
        fid = 1;
    else
        logging = false;
    end

    filenameBusData_Area = strcat(systemDataFolder, "area", num2str(Area), filesep, 'powerdata.csv');
    busData_Area = readmatrix(filenameBusData_Area);
    busDataTable_Area = array2table(busData_Area, 'VariableNames', {'bus', 'P_L', 'Q_L', 'Q_C', 'P_der', 'busType'});
    
    N_Area = size(busData_Area, 1);     % total number of nodes
    
    filenameBranchData_Area = strcat(systemDataFolder, "area", num2str(Area), '/linedata.csv');
    branchData_Area = readmatrix(filenameBranchData_Area);
    branchDataTable_Area = array2table(branchData_Area, 'VariableNames', {'fb', 'tb', 'R', 'X'});

    filenameActualBusNums_Area = strcat(systemDataFolder, "area", num2str(Area), '/powerdataActualBusNums.csv');
    opts = detectImportOptions(filenameActualBusNums_Area);
    opts.DataLines = 2;
    opts.VariableNamesLine = 1;
    busDataActualBusNumsTable_Area = readtable(filenameActualBusNums_Area, opts);
    
    if displayTables
        display(busDataTable_Area)
        display(branchDataTable_Area)
        display(busDataActualBusNumsTable_Area)
    end
    
    writetable(sortrows(busDataTable_Area, "bus"))
    % Graph Formation
    fb_Area = branchDataTable_Area.fb;
    tb_Area = branchDataTable_Area.tb;
    graph_Area = graph(fb_Area, tb_Area);
    edgeMatrix_Area = [fb_Area, tb_Area];
    
% Optional: Plot Graphs which highlight the relationships between different Areas.

    if macroItr == 1 && t == 1
        plotGraphs(t, macroItr, Area, ...
        N_Area, ...
        busDataActualBusNumsTable_Area, graph_Area, ...
        isRoot_Area, numAreas, CB_FullTable, ...
        numChildAreas_Area, 'verbose', verbose, 'logging', logging, 'displayNetworkGraphs', displayNetworkGraphs, 'displayActualBusNumbersInGraphs', displayActualBusNumbersInGraphs, 'savePlots', savePlots);
    end
    
    % Line Data
    Z_B = kV_B^2*1000/kVA_B;                           % base Z
    
    R_Area = branchDataTable_Area.R/Z_B;                     % R of edge
    X_Area = branchDataTable_Area.X/Z_B;                     % X of edge
    
    P_L_Area = busDataTable_Area.P_L.*load_mult/kVA_B;           % Rated Pload of node
    Q_L_Area = busDataTable_Area.Q_L.*load_mult/kVA_B;           % Rated Qload of node
    Q_C_Area = busDataTable_Area.Q_C/kVA_B;                           % Rated Q of Capacitor
    
    P_der_Area = busDataTable_Area.P_der*gen_mult/kVA_B;          % Rated  DER active power
    S_der_Area = 1.2*busDataTable_Area.P_der/kVA_B; 
    
    busDataTable_pu_Area = busDataTable_Area;
    busDataTable_pu_Area.P_L = P_L_Area;
    busDataTable_pu_Area.Q_L = Q_L_Area;
    busDataTable_pu_Area.Q_C = Q_C_Area;
    busDataTable_pu_Area.P_der = P_der_Area;
    busDataTable_pu_Area.S_der = S_der_Area;

end