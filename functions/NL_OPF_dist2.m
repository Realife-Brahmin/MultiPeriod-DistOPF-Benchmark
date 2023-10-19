function [x, sysInfo, simInfo, ...
     time_dist] = ...
    ...
    NL_OPF_dist2(sysInfo, simInfo, areaInfo, v_parent_Area_1toT, S_connection_Area_1toT,  ...
    lambdaVals, pvCoeffVals, time_dist, varargin)
    
 % Default values for optional arguments
    verbose = false;
    logging_Aeq_beq = false;
    logging = false;
    CVR = [0; 0];
    V_max = 1.05;
    V_min = 0.95;
    Qref_DER = 0.00;
    Vref_DER = 1.00;
    delta_t = 0.25;
    etta_D = 0.95;
    etta_C = 0.95;
    chargeToPowerRatio = 4;
    soc_min = 0.30;
    soc_max = 0.95;
    alpha = 1e-3;
    gamma = 1e0;
    profiling = false;
    saveSCDPlots = false;
    displayTables = false;

    saveToFile = false;
    Area = areaInfo.Area;
    strArea = convert2doubleDigits(Area);
    fileExtension = ".txt";
    systemName = "ieee123";
    saveLocationName = "logfiles/";
    fileOpenedFlag_Aeq_beq = false;

    % Process optional arguments
    numArgs = numel(varargin);

    if mod(numArgs, 2) ~= 0
        error('Optional arguments must be specified as name-value pairs.');
    end
    
    validArgs = ["verbose", "logging", "logging_Aeq_beq", "systemName", ...
        "CVR", "V_max", "V_min", "saveToFile", "saveLocation", "Qref_DER", ...
        "Vref_DER", "fileExtension", "delta_t", "etta_C", "etta_D", ...
        "chargeToPowerRatio", "soc_min", "soc_max", "profiling", ...
        "saveSCDPlots", "displayTables"];
    
    for kwarg_num = 1:2:numArgs
        argName = varargin{kwarg_num};
        argValue = varargin{kwarg_num+1};
        
        if ~ischar(argName) || ~any(argName == validArgs)
            error('Invalid optional argument name.');
        end
        
        switch argName
            case "verbose"
                verbose = argValue;
            case "logging"
                logging = argValue;
            case "logging_Aeq_beq"
                logging_Aeq_beq = argValue;
            case "systemName"
                systemName = argValue;
            case "CVR"
                CVR = argValue;
            case "V_max"
                V_max = argValue;
            case "V_min"
                V_min = argValue;
            case 'saveToFile'
                saveToFile = argValue;
            case 'saveLocation'
                saveLocationName = argValue;
            case 'Vref_DER'
                Vref_DER = argValue;
            case 'Qref_DER'
                Qref_DER = argValue;
            case 'delta_t'
                delta_t = argValue;
            case 'etta_C'
                etta_C = argValue;
            case 'etta_D'
                etta_D = argValue;
            case "chargeToPowerRatio"
                chargeToPowerRatio = argValue;
            case "soc_min"
                soc_min = argValue;
            case "soc_max"
                soc_max = argValue;
            case "alpha"
                alpha = argValue;
            case "gamma"
                gamma = argValue;
            case "profiling"
                profiling = argValue;
            case "saveSCDPlots"
                saveSCDPlots = argValue;
            case "displayTables"
                displayTables = argValue;
        end
    end
    
    numAreas = sysInfo.numAreas;

    saveLocationFilename = strcat(saveLocationName , systemName, "/numAreas_", num2str(numAreas), "/optimizationLogs", fileExtension);
    saveLocationFilename_Aeq_beq = strcat(saveLocationName , systemName, "/numAreas_", num2str(numAreas), "/Aeq_beq_area", strArea, fileExtension);
    
    macroItr = simInfo.macroItr; % completed macro-iterations, starts at 0

    if macroItr ~= 0
        logging_Aeq_beq = false;
    end

    if logging_Aeq_beq && saveToFile && macroItr == 0 && Area == 2
        fileOpenedFlag_Aeq_beq = true;
        fid_Aeq_beq = fopen(saveLocationFilename_Aeq_beq, 'w');  % Open file for writing
    else
        logging_Aeq_beq = false;
        fid_Aeq_beq = 1;
    end
    
    if logging && verbose
        error("Kindly specify ONLY one of the following arguments as true: verbose and logging.")
    elseif logging && ~verbose
        fileOpenedFlag = true;
        fid = fopen(saveLocationFilename, 'w');
    elseif ~logging
        logging = verbose;
        fid = 1;
    end
    
    logging_Aeq_beq = false;
    
    isRoot_Area = sysInfo.isRoot(Area);
    numAreas = sysInfo.numAreas;
    systemName = sysInfo.systemName;
    CB_FullTable = sysInfo.CBTable;
    numChildAreas_Area = sysInfo.numChildAreas(Area);

    areaInfo = extractAreaInfo(areaInfo, simInfo, isRoot_Area, systemName, numAreas, ...
    CB_FullTable, numChildAreas_Area, 'verbose', verbose, 'logging', logging, 'displayNetworkGraphs', false, 'displayTables', displayTables);
    
    areaInfo = exchangeCompVars(areaInfo, S_connection_Area_1toT);
    
     myfprintf(logging_Aeq_beq, fid_Aeq_beq, "**********" + ...
        "Constructing Aeq and beq for Area %d.\n" + ...
        "***********\n", Area); 

    CVR_P = CVR(1);
    CVR_Q = CVR(2);
    
    [Aeq, beq, lb, ub, x0, areaInfo] = LinEqualities(areaInfo, simInfo, lambdaVals, pvCoeffVals, v_parent_Area_1toT);

    % plotSparsity(Aeq, beq);
    
    t3Start = tic;
    
    microItrMax = simInfo.alg.microItrMax;
    tolfun = simInfo.alg.tolfun;
    stepTol = simInfo.alg.stepTol;
    constraintTol = simInfo.alg.constraintTol;
    optimalityTol = simInfo.alg.optimalityTol;
    % printouts = 'iter-detailed';
    printouts = 'off';
    % options = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxIterations', microItrMax, 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp');
    options = optimoptions('fmincon', 'Display', printouts, 'MaxIterations', microItrMax, 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp', ...
        'FunctionTolerance', tolfun, ...
        'StepTolerance', stepTol, ...             % Equivalent to 10 watts
    'ConstraintTolerance', constraintTol, ...       % For line flows, voltages, etc.
    'OptimalityTolerance', optimalityTol);
    % options = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxIterations', microItrMax, 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp', 'PlotFcn', @optimplotfval);
    % options = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxIterations', 200, 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp', 'PlotFcn', @optimplotfval);
    
    T = simInfo.T;

    if T >= 7
        profiling = true;
        profile on
    end
    objectiveFuns = {"func_PLoss", "func_SCD", "func_netChangeInSOC"};
    % objectiveFuns = {"func_PLoss", "func_SCD"};
    % objectiveFuns = {"func_PLoss", "func_netChangeInSOC"};
    [x, fval, ~, output] = fmincon(@(x)objfun(x, areaInfo, T, 'objectiveFuns', objectiveFuns, 'alpha', alpha, 'gamma', gamma), ...
    x0, [], [], Aeq, beq, lb, ub, ...
    @(x)NonLinEqualities(x, areaInfo, T, "verbose", false, "saveToFile", false), ...
    options);

    iterations_taken = output.iterations;
    
    folderName = strcat("processedData", filesep, sysInfo.systemName, filesep, "numAreas_", num2str(numAreas), filesep, "area", num2str(Area));
    if ~exist(folderName, 'dir')
        mkdir(folderName)
    end
    prefixName = strcat(folderName, filesep, "Horizon_", num2str(T), "_macroItr_", num2str(macroItr+1));
    areaSolutionName_x = strcat(prefixName, "_optimalSolutions.csv");
    areaSolutionName_fval = strcat(prefixName, "_", getenv('COMPUTERNAME'), "_optimalObjectiveFunctionValue.txt");
    
    writematrix(x, areaSolutionName_x);

    if profiling
        profile viewer;
    end
    lineLosses = objfun(x, areaInfo, T, 'objectiveFuns', {"func_PLoss"});
    scd = objfun(x, areaInfo, T, 'objectiveFuns', {"func_SCD"});
    changeInSOC = objfun(x, areaInfo, T, 'objectiveFuns', {"func_netChangeInSOC"});
    
    t3 = toc(t3Start);

    % filename = sprintf('%s', areaSolutionName_fval);

    % Open the file for writing
    fid = fopen(areaSolutionName_fval, 'w');
    fileOpenedFlag = true;

    % Check if the file opened successfully
    if fid == -1
        error('Failed to open the file for writing.');
    end

    [nLinEqnsT, nNonLinEqnsT, nVarsT, areaInfo] = getProblemSize(areaInfo, T);
    
    sysInfo.Area{Area} = areaInfo;

    myfprintf(true, fid, "Machine this simulation was solved on: %s\n", getenv('COMPUTERNAME'));
    myfprintf(true, fid, "Optimization for Area %d for %d time periods took %d [s] and %d iterations.\n", Area, T,  t3, iterations_taken);
    myfprintf(true, fid, "where a micro-iteration limit of %d and a minimum improvement of %d [kW] was imposed.\n", microItrMax, tolfun*1e3);
    myfprintf(true, fid, "Average time per iteration: %d [s]\n", t3/iterations_taken);
    myfprintf(true, fid, "Number of Linear Equations: %d\n", nLinEqnsT);
    myfprintf(true, fid, "Number of Nonlinear Equalities: %d\n", nNonLinEqnsT);
    myfprintf(true, fid, "Number of Optimization Variables: %d\n", nVarsT);
    myfprintf(true, fid, "Total Objective Function Value for Area %d for %d time periods = %d [kW]\n", Area, T, fval*1000);
    myfprintf(true, fid, "Total Real Power Line Losses for Area %d for %d time periods = %d [kW]\n", Area, T, lineLosses*1000);
    myfprintf(true, fid, "Total Battery Power losses for Area %d for %d time periods = %d [kW]\n", Area, T, scd*1000/alpha);
    myfprintf(true, fid, "Average SOC Level constraint violation for Area %d for %d time periods = %d [kWh]\n", Area, T, sqrt( changeInSOC*(1000^2)/(gamma * areaInfo.nBatt_Area)));
    myfprintf(true, fid, "where alpha = %d and gamma = %d\n", alpha, gamma);
   
    saveSCDPlots = ~macroItr && saveSCDPlots;

    if saveSCDPlots
        checkForSCD(false, sysInfo, simInfo, areaInfo, T, x, 'savePlots', saveSCDPlots);
    end
    
    time_dist(macroItr+1, Area) = t3;


    if fileOpenedFlag
        fclose(fid);
    end

end