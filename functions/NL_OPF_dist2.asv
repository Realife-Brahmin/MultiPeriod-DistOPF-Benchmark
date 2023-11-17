function [x, sysInfo, simInfo, ...
     time_dist] = ...
    ...
    NL_OPF_dist2(sysInfo, simInfo, areaInfo, v_parAr_1toT, S_chArs_1toT,  ...
    time_dist, varargin)
    
    noBatteries = simInfo.noBatteries;
 % Default values for optional arguments
    verbose = false;
    logging_Aeq_beq = false;
    logging = false;
    CVR = [0; 0];
    % % V_max = 1.05;
    % V_max = simInfo.V_max;
    % % V_min = 0.95;
    % V_min = simInfo.V_min;
    Qref_DER = 0.00;
    Vref_DER = 1.00;
    % % delta_t = 0.25;
    % delta_t = simInfo.delta_t;
    % % etta_C = 0.95;
    % etta_C = simInfo.etta_C;
    % % etta_D = 0.95;
    % etta_D = simInfo.etta_D;
    % % chargeToPowerRatio = 4;
    % chargeToPowerRatio = simInfo.chargeToPowerRatio;
    % % soc_min = 0.30;
    % soc_min = simInfo.soc_min;
    % % soc_max = 0.95;
    % soc_max = simInfo.soc_max;
    % alpha = 1e-3;
    alpha = simInfo.alpha;
    % gamma = 1e0;
    % gamma = 1e-1;
    gamma = simInfo.gamma;
    profiling = false;
    saveSCDPlots = false;
    displayTables = false;
    plotLinDS = false;

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
        "CVR", "saveToFile", "saveLocation", "Qref_DER", ...
        "Vref_DER", "fileExtension", ...
        "profiling", ...
        "saveSCDPlots", "displayTables", "plotLinDS"];
    
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
            case 'saveToFile'
                saveToFile = argValue;
            case 'saveLocation'
                saveLocationName = argValue;
            case 'Vref_DER'
                Vref_DER = argValue;
            case 'Qref_DER'
                Qref_DER = argValue;
            case "profiling"
                profiling = argValue;
            case "saveSCDPlots"
                saveSCDPlots = argValue;
            case "displayTables"
                displayTables = argValue;
            case "plotLinDS"
                plotLinDS = argValue;
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
    
    areaInfo = exchangeCompVars(areaInfo, S_chArs_1toT);
    
    myfprintf(logging_Aeq_beq, fid_Aeq_beq, "**********" + ...
        "Constructing Aeq and beq for Area %d.\n" + ...
        "***********\n", Area); 

    % CVR_P = CVR(1);
    % CVR_Q = CVR(2);
    
    [Aeq, beq, lb, ub, x0, areaInfo] = LinEqualities(areaInfo, simInfo, v_parAr_1toT);
    
    % plotLinDS = true;

    if plotLinDS && macroItr == 0
        plotSparsity(Aeq, beq);
    end
    
    ext = ".csv";
    if ~noBatteries
        battstring = "withBatteries";
    else
        battstring = "withoutBatteries";
    end
    saveLocationFolderName = strcat("processedData", filesep , systemName, filesep, "numAreas_", num2str(numAreas), filesep, "Area", num2str(Area));
    if ~exist("saveLocationFolderName", 'dir')
        mkdir(saveLocationFolderName)
    end
    filenameAeq = strcat(saveLocationFolderName, filesep, "Aeq_B_T_", num2str(simInfo.T), "_", battstring, "_macroItr_", num2str(1+macroItr), ext);
    % Saving Aeq to a CSV file
    writematrix(Aeq, filenameAeq);
    filenamebeq = strcat(saveLocationFolderName, filesep, "beq_B_T_", num2str(simInfo.T), "_", battstring,  "_macroItr_", num2str(1+macroItr), ext);

    % Saving beq to a CSV file
    writematrix(beq, filenamebeq);

    t3Start = tic;
    
    microItrMax = simInfo.alg.microItrMax;
    tolfun = simInfo.alg.tolfun;
    stepTol = simInfo.alg.stepTol;
    constraintTol = simInfo.alg.constraintTol;
    optimalityTol = simInfo.alg.optimalityTol;
    displayIterations = 'off';
    
    options = optimoptions('fmincon', 'Display', displayIterations, 'MaxIterations', microItrMax, 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp', ...
        'FunctionTolerance', tolfun, ...
        'StepTolerance', stepTol, ...             % Equivalent to 10 watts
    'ConstraintTolerance', constraintTol, ...       % For line flows, voltages, etc.
    'OptimalityTolerance', optimalityTol);
    
    T = simInfo.T;

    if T >= 7
        profiling = true;
        profile on
    end
    
    objFunction = simInfo.objFunction;
    if strcmp(objFunction, "loss_min")
        mainObjFunc = "func_PLoss";
    elseif strcmp(objFunction, "gen_cost")
        mainObjFunc = "func_PSubsCost";
    else
        error("Unknown objective function");
    end

    if ~noBatteries
        objectiveFuns = {mainObjFunc, "func_SCD", "func_netChangeInSOC"};
        % objectiveFuns = {mainObjFunc, "func_SCD"};
        % objectiveFuns = {mainObjFunc, "func_netChangeInSOC"};
    else
        objectiveFuns = {mainObjFunc};
    end
    
    % lb
    % ub
    [x, fval, ~, output] = fmincon(@(x)objfun(x, simInfo, sysInfo, areaInfo, T, 'objectiveFuns', objectiveFuns), ...
    x0, [], [], Aeq, beq, lb, ub, ...
    @(x)NonLinEqualities(x, simInfo, areaInfo, T, "verbose", false, "saveToFile", false), ...
    options);
    
    iterations_taken = output.iterations;
    
    folderName = strcat("processedData", filesep, sysInfo.systemName, filesep, "numAreas_", num2str(numAreas), filesep, "area", num2str(Area));
    if ~exist(folderName, 'dir')
        mkdir(folderName)
    end
    prefixName = strcat(folderName, filesep, "Horizon_", num2str(T), "_macroItr_", num2str(macroItr+1));
    areaSolutionName_x = strcat(prefixName, "_optimalSolutions.csv");
    if ~noBatteries
        battstring = "withBatteries";
    else
        battstring = "withoutBatteries";
    end
    areaSolutionName_fval = strcat(prefixName, "_", getenv('COMPUTERNAME'), "_optimalObjectiveFunctionValue_", battstring, ".txt");
    
    writematrix(x, areaSolutionName_x);

    if profiling
        profile viewer;
    end

    [lineLosses, ~] = objfun(x, simInfo, sysInfo, areaInfo, T, 'objectiveFuns', {"func_PLoss"});
    [~, lineLosses_1toT] = objfun(x, simInfo, sysInfo, areaInfo, T, 'objectiveFuns', {"func_PLoss_1toT"});
    [genPower, ~] = objfun(x, simInfo, sysInfo, areaInfo, T, 'objectiveFuns', {"func_PSubs"});
    [~, genPower_1toT] = objfun(x, simInfo, sysInfo, areaInfo, T, 'objectiveFuns', {"func_PSubs_1toT"});

    [genCost, ~] = objfun(x, simInfo, sysInfo, areaInfo, T, 'objectiveFuns', {"func_PSubsCost"});
    [~, genCost_1toT] = objfun(x, simInfo, sysInfo, areaInfo, T, 'objectiveFuns', {"func_PSubsCost_1toT"});

    if ~noBatteries
        [scd, ~] = objfun(x, simInfo, sysInfo, areaInfo, T, 'objectiveFuns', {"func_SCD"});
        [changeInSOC, ~] = objfun(x, simInfo, sysInfo, areaInfo, T, 'objectiveFuns', {"func_netChangeInSOC"});
    end

    t3 = toc(t3Start);

    % Open the file for writing
    fid = fopen(areaSolutionName_fval, 'w');
    fileOpenedFlag = true;

    % Check if the file opened successfully
    if fid == -1
        error('Failed to open the file for writing.');
    end

    [nLinEqnsT, nNonLinEqnsT, nVarsT, areaInfo] = getProblemSize(areaInfo, T);
    kVA_B = sysInfo.kVA_B;
    kV_B = sysInfo.kV_B;

    myfprintf(true, fid, "Machine this simulation was solved on: %s\n", getenv('COMPUTERNAME'));
    myfprintf(true, fid, "Optimization for Area %d for %d time periods took %d [s] and %d iterations.\n", Area, T,  t3, iterations_taken);
    myfprintf(true, fid, "where a micro-iteration limit of %d and a minimum improvement of %d [kW] was imposed.\n", microItrMax, tolfun*1e3);
    myfprintf(true, fid, "Average time per iteration: %d [s]\n", t3/iterations_taken);
    myfprintf(true, fid, "Number of Linear Equations: %d\n", nLinEqnsT);
    myfprintf(true, fid, "Number of Nonlinear Equalities: %d\n", nNonLinEqnsT);
    myfprintf(true, fid, "Number of Optimization Variables: %d\n", nVarsT);
    myfprintf(true, fid, "Total Objective Function Value for Area %d for %d time periods = %d \n", Area, T, fval);
    myfprintf(true, fid, "Total Real Power Line Losses for Area %d for %d time periods = %d [kW]\n", Area, T, lineLosses*kVA_B);
    if ~noBatteries
        myfprintf(true, fid, "Total Battery Power losses for Area %d for %d time periods = %d [kW]\n", Area, T, scd*kVA_B/alpha);
        myfprintf(true, fid, "Average SOC Level constraint violation for Area %d for %d time periods = %d [kWh]\n", Area, T, sqrt( changeInSOC*(kVA_B^2)/(gamma * areaInfo.nBatt_Area)));
        myfprintf(true, fid, "where alpha = %d and gamma = %d\n", alpha, gamma);
    end
    
    saveSCDPlots = ~macroItr && saveSCDPlots;
    % saveSCDPlots = false;
    % saveSCDPlots = true;
    % keyboard;
    if ~noBatteries && saveSCDPlots
        checkForSCD(sysInfo, simInfo, areaInfo, T, x, 'savePlots', true);
    end
    
    time_dist(macroItr+1, Area) = t3;
    
    N_Area = areaInfo.N_Area;
    m_Area = areaInfo.m_Area;
    nDER_Area = areaInfo.nDER_Area;
    nBatt_Area = areaInfo.nBatt_Area;
    
    xVals_Area = x;
    P_Area_1toT = reshape(xVals_Area(areaInfo.indices_Pij), m_Area, T); %m_Areax1
    Q_Area_1toT = reshape(xVals_Area(areaInfo.indices_Qij), m_Area, T); %m_Areax1
    S_Area_1toT = complex(P_Area_1toT, Q_Area_1toT); %m_Areax1

    l_Area_1toT = reshape(xVals_Area(areaInfo.indices_lij), m_Area, T);
    vAll_Area_1toT = reshape(xVals_Area(areaInfo.indices_vAllj), N_Area, T); %N_Areax1
    
    qD_Area_1toT = reshape(xVals_Area(areaInfo.indices_qDj), nDER_Area, T);
    
    if ~noBatteries
        B_Area_1toT = reshape(xVals_Area(areaInfo.indices_Bj), nBatt_Area, T);
        Pc_Area_1toT = reshape(xVals_Area(areaInfo.indices_Pcj), nBatt_Area, T);
        Pd_Area_1toT = reshape(xVals_Area(areaInfo.indices_Pdj), nBatt_Area, T);
        qB_Area_1toT = reshape(xVals_Area(areaInfo.indices_qBj), nBatt_Area, T);

        areaInfo.B_Area_1toT = B_Area_1toT;
        areaInfo.Pc_Area_1toT = Pc_Area_1toT;
        areaInfo.Pd_Area_1toT = Pd_Area_1toT;
        areaInfo.qB_Area_1toT = qB_Area_1toT;

        areaInfo.scd = scd;
        areaInfo.changeInSOC = changeInSOC;
        areaInfo.genCost = genCost;
    else

        areaInfo.B_Area_1toT = "NA";
        areaInfo.Pc_Area_1toT = "NA";
        areaInfo.Pd_Area_1toT = "NA";
        areaInfo.qB_Area_1toT = "NA";

        areaInfo.scd = "NA";
        areaInfo.changeInSOC = "NA";
    end

    areaInfo.P_Area_1toT = P_Area_1toT;
    areaInfo.Q_Area_1toT = Q_Area_1toT;
    areaInfo.S_Area_1toT = S_Area_1toT;

    areaInfo.l_Area_1toT = l_Area_1toT;
    areaInfo.vAll_Area_1toT = vAll_Area_1toT;

    areaInfo.qD_Area_1toT = qD_Area_1toT;
    
    
    PLoss_allT = lineLosses;
    PLoss_1toT = lineLosses_1toT;
    areaInfo.PLoss_allT = PLoss_allT;
    areaInfo.PLoss_1toT = PLoss_1toT;
    PSubs_allT = genPower;
    PSubs_1toT = genPower_1toT;
    areaInfo.PSubs_allT = PSubs_allT;
    areaInfo.PSubs_1toT = PSubs_1toT;
    PSubsCost_allT = genCost;
    PSubsCost_1toT = genCost_1toT;
    areaInfo.PSubsCost_allT = PSubsCost_allT;
    areaInfo.PSubsCost_1toT = PSubsCost_1toT;
    areaInfo.fval = fval;
    areaInfo.xvals = xVals_Area;
    sysInfo.Area{Area} = areaInfo;

    if fileOpenedFlag
        fclose(fid);
    end

end