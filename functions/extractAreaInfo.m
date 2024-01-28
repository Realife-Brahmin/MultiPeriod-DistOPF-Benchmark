function areaInfo = extractAreaInfo(areaInfo, sysInfo, simInfo, isRoot_Area, systemName, numAreas, CB_FullTable, numChildAreas_Area, varargin)

 % Default values for optional arguments
    verbose = false;
    logging = false;
    displayTables = false;
    displayNetworkGraphs = true;
    displayActualBusNumbersInGraphs = false;
    saveLocationName = "logfiles/";
    fileExtension = ".txt";
    savePlots = true;
    kV_B = sysInfo.kV_B;
    kVA_B = sysInfo.kVA_B;
    DER_percent = simInfo.DER_percent;
    S_to_P_ratio_PV = simInfo.S_to_P_ratio_PV;
    Batt_percent = simInfo.Batt_percent;

    load_mult = 1;
    gen_mult = 1;
    systemDataFolder = strcat("rawData", filesep, systemName, filesep, "numAreas_", num2str(numAreas), filesep);

    % Process optional arguments
    numArgs = numel(varargin);
    if mod(numArgs, 2) ~= 0
        error('Optional arguments must be specified as name-value pairs.');
    end
    
    validArgs = ["verbose", "logging", "displayTables", "displayNetworkGraphs", ...
        "displayActualNumbersInGraphs", "savePlots", "load_mult", "gen_mult"];
    
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
    
    T = simInfo.T;
    macroItr = simInfo.macroItr; % completed macro-iterations, starts at 0
    lambdaVals = simInfo.lambdaVals;
    pvCoeffVals = simInfo.pvCoeffVals;

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
    
    filenameBusData_Area = strcat(systemDataFolder, "area", num2str(Area), filesep, 'powerdata.csv');
    busData_Area = readmatrix(filenameBusData_Area);
    busDataTable_Area = array2table(busData_Area, 'VariableNames', {'bus', 'P_L', 'Q_L', 'Q_C', 'P_der', 'busType'});
    % % display(busDataTable_Area)
    Area = areaInfo.Area;
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
    
    writetable(sortrows(busDataTable_Area, "bus"), "sortedBusData.csv")
    % Graph Formation
    fb_Area = branchDataTable_Area.fb;
    tb_Area = branchDataTable_Area.tb;
    graph_Area = graph(fb_Area, tb_Area);
    
    % Line Data
    MVA_B = kVA_B * 1e-3;

    Z_B = kV_B^2/MVA_B;

    R_Area = branchDataTable_Area.R/Z_B;                     % R of edge
    X_Area = branchDataTable_Area.X/Z_B;                     % X of edge
    
    P_L_Area = busDataTable_Area.P_L.*load_mult/kVA_B;           % Rated Pload of node
    Q_L_Area = busDataTable_Area.Q_L.*load_mult/kVA_B;           % Rated Qload of node
    Q_C_Area = busDataTable_Area.Q_C/kVA_B;                           % Rated Q of Capacitor
    
    P_der_Area0 = busDataTable_Area.P_der*gen_mult/kVA_B;          % Rated  DER active power, all DERs from original CSV file
    P_batt_Area0 = busDataTable_Area.P_der*gen_mult/kVA_B; % Rated Battery Active Power
    P_der_Area = select_percentage_of_nz_elements(P_der_Area0, DER_percent);          % Rated  DER active power, only DER_percent elements selected
    P_batt_Area = select_percentage_of_nz_elements(P_batt_Area0, Batt_percent); % Rated Battery Real Power, only Batt_percent elements selected

    S_der_Area0 = S_to_P_ratio_PV*busDataTable_Area.P_der/kVA_B; 
    S_der_Area = select_percentage_of_nz_elements(S_der_Area0, DER_percent);
    S_batt_Area0 = S_to_P_ratio_PV*busDataTable_Area.P_der/kVA_B;
    S_batt_Area = select_percentage_of_nz_elements(S_batt_Area0, Batt_percent);

    P_L_Area_1toT = repmat(P_L_Area, 1, T).*lambdaVals;
    Q_L_Area_1toT = repmat(Q_L_Area, 1, T).*lambdaVals;

    P_der_Area0_1toT = repmat(P_der_Area0, 1, T).*pvCoeffVals;
    P_der_Area_1toT = repmat(P_der_Area, 1, T).*pvCoeffVals;

    busDataTable_pu_Area = tableToStructColumnwise(busDataTable_Area);
    busDataTable_pu_Area.P_L_1toT = P_L_Area_1toT;
    busDataTable_pu_Area.Q_L_1toT = Q_L_Area_1toT;
    busDataTable_pu_Area.Q_C = Q_C_Area;
    busDataTable_pu_Area.P_der0_1toT = P_der_Area0_1toT;
    busDataTable_pu_Area.P_batt0 = P_batt_Area0;
    busDataTable_pu_Area.P_der_1toT = P_der_Area_1toT;
    busDataTable_pu_Area.P_batt = P_batt_Area;
    busDataTable_pu_Area.S_der0 = S_der_Area0;
    busDataTable_pu_Area.S_batt0 = S_batt_Area0;

    busDataTable_pu_Area.S_der = S_der_Area;
    busDataTable_pu_Area.S_batt = S_batt_Area;

    areaInfo = getAreaParameters(simInfo, Area, busDataTable_pu_Area, branchDataTable_Area, R_Area, X_Area);

    % Optional: Plot Graphs which highlight the relationships between different Areas.
    
    if macroItr == 0
        plotGraphs(areaInfo, simInfo, ...
        busDataActualBusNumsTable_Area, graph_Area, ...
        isRoot_Area, numAreas, CB_FullTable, ...
        numChildAreas_Area, 'verbose', verbose, 'logging', logging, ...
        'displayNetworkGraphs', displayNetworkGraphs, ...
        'displayActualBusNumbersInGraphs', displayActualBusNumbersInGraphs, ...
        'savePlots', savePlots);
    end
    

end