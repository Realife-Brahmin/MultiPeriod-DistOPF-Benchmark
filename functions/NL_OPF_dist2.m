function [x, B0Vals_pu_Area, ...
    macroIterationPLosses, macroIterationQLosses, ...
    macroIterationPSaves, macroItr, time_dist, ...
    N_Area, m_Area, nDER_Area, nBatt_Area, ...
    busDataTable_Area, branchDataTable_Area] = ...
    ...
    NL_OPF_dist2(v_parent_Area, S_connection_Area, B0Vals_pu_Area, ...
    lambdaVals, pvCoeffVals, ...
    Area, isLeaf_Area, isRoot_Area, numChildAreas_Area, numAreas, ...
    macroIterationPLosses, macroIterationQLosses, macroIterationPSaves, ...
    macroItr, time_dist, t, T, ...
    CB_FullTable, varargin)
    
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

    saveToFile = false;
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
        "chargeToPowerRatio", "soc_min", "soc_max"];
    
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
        end
    end
    
    saveLocationFilename = strcat(saveLocationName , systemName, "/numAreas_", num2str(numAreas), "/optimizationLogs", fileExtension);
    saveLocationFilename_Aeq_beq = strcat(saveLocationName , systemName, "/numAreas_", num2str(numAreas), "/Aeq_beq_area", strArea, fileExtension);

    if macroItr ~= 1
        logging_Aeq_beq = false;
    end

    if logging_Aeq_beq && saveToFile && macroItr == 1 && Area == 2
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
        if t == 1
            fid = fopen(saveLocationFilename, 'w');
        else
            fid = fopen(saveLocationFilename, 'a');
        end
    elseif ~logging
        logging = verbose;
        fid = 1;
    end
    
    logging_Aeq_beq = false;

    [busDataTable_Area, branchDataTable_Area, edgeMatrix_Area, R_Area, X_Area] ...
        = extractAreaElectricalParameters(Area, t, macroItr, isRoot_Area, systemName, numAreas, CB_FullTable, numChildAreas_Area, 'verbose', verbose, 'logging', logging, 'displayNetworkGraphs', false);
    
    areaInfo = getAreaParameters(Area, busDataTable_Area, branchDataTable_Area, R_Area, X_Area);
    areaInfo = exchangeCompVars(areaInfo, S_connection_Area);

    
     myfprintf(logging_Aeq_beq, fid_Aeq_beq, "**********" + ...
        "Constructing Aeq and beq for Area %d.\n" + ...
        "***********\n", Area); 

    CVR_P = CVR(1);
    CVR_Q = CVR(2);
    
    [Aeq, beq, lb, ub, x0, areaInfo] = LinEqualities(areaInfo, T, lambdaVals, pvCoeffVals, v_parent_Area);

    % plotSparsity(Aeq, beq);
    
    t3Start = tic;

    function stop = outfun(x, optimValues, state)
        stop = false;
        disp(['Current X: ' num2str(x)]);
        disp(['Current function value: ' num2str(optimValues.fval)]);
    end
    
    % itermax = 50;
    itermax = 20;
    % options = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxIterations', 20, 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp');
    options = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxIterations', itermax, 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp', 'PlotFcn', @optimplotfval);
    % options = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxIterations', 200, 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp', 'PlotFcn', @optimplotfval);

    profile on
    % objectiveFuns = {"func_PLoss", "func_SCD", "func_netChangeInSOC"};
    objectiveFuns = {"func_PLoss", "func_SCD"};
    % objectiveFuns = {"func_PLoss", "func_netChangeInSOC"};
    try
        [x, fval] = fmincon( @(x)objfun(x, areaInfo, T, 'objectiveFuns', objectiveFuns), ...
            x0, [], [], Aeq, beq, lb, ub, ...
            @(x)NonLinEqualities(x, areaInfo, T, "verbose", false, "saveToFile", false), ...
            options);

        profile viewer;
        keyboard;

        lineLosses = objfun(x, areaInfo, T, 'objectiveFuns', {"func_PLoss"});
        scd = objfun(x, areaInfo, T, 'objectiveFuns', {"func_SCD"});
        changeInSOC = objfun(x, areaInfo, T, 'objectiveFuns', {"func_netChangeInSOC"});

        myfprintf(true, "Real Power Line Losses for Area %d for %d time periods = %d [kW]\n", Area, T, lineLosses*1000);
        myfprintf(true, "SCD Constraint violation for Area %d for %d time periods = %d [kW]\n", Area, T, scd*1000);
        myfprintf(true, "SOC Level constraint violation for Area %d for %d time period = %d [kWh]\n", Area, T, changeInSOC*1000);
        checkForSCD(areaInfo, T, x)
        keyboard;
    catch ME
        % Save the current workspace to a file for debugging
        save('errorWorkspace.mat');
        rethrow(ME);
    end
    % macroIterationPLoss = fval;
    macroIterationQLoss = objfun(x, areaInfo, T, 'objectiveFuns', {"func_QLoss"});

    
    t3 = toc(t3Start);
    
    time_dist(macroItr, Area) = t3;
    
    % % Result
    % P_Area = x(indices_Pij_T); %m_Areax1
    % Q_Area = x(indices_Q); %m_Areax1
    % % S_Area = complex(P_Area, Q_Area); %m_Areax1
    % % l_Area = x(indices_l);
    % % v_Area = x(indices_vAll); %N_Areax1
    % qD_Area = x(indices_qD);
    % % v_Area(1) = v_parent_Area;
    % B_Area = x(indices_B);
    % Pd_Area = x(indices_Pd);
    % Pc_Area = x(indices_Pc);
    % qB_Area = x(indices_qB);
    % 
    % qD_AllBuses = zeros(N_Area, 1);
    % 
    % for der_num = 1 : nDER_Area
    %     busNum = busesWithDERs_Area(der_num);
    %     qD_AllBuses(busNum) = qD_Area(der_num);
    % end
    % 
    % [B_AllBuses, Pd_AllBuses, Pc_AllBuses, qB_AllBuses] = deal(zeros(N_Area, 1));
    % 
    % for batt_num = 1 : nBatt_Area
    %     busNum = busesWithBatts_Area(batt_num); 
    %     B_AllBuses(busNum) = B_Area(batt_num);
    %     Pc_AllBuses(busNum) = Pc_Area(batt_num);
    %     Pd_AllBuses(busNum) = Pd_Area(batt_num);
    %     qB_AllBuses(busNum) = qB_Area(batt_num);
    % end
    
    % P_inFlowArea = P_Area(1);
    % P_der_Total = sum(P_der_Area);
    % Pd_Total = sum(Pd_Area);
    % Pc_Total = sum(Pc_Area);
    % PLoad_Total = sum(P_L_Area);
    % 
    % PLoss = P_inFlowArea + P_der_Total + Pd_Total - PLoad_Total - Pc_Total;
    % percentageSavings = 100* (Pd_Total - Pc_Total) / (P_inFlowArea + P_der_Total + Pd_Total - Pc_Total);
    
    % Q_inFlowArea = Q_Area(1);
    % qD_Total = sum(qD_Area);
    % qB_Total = sum(qB_Area);
    % QLoad_Total = sum(Q_L_Area);
    % 
    % QLoss = Q_inFlowArea + qD_Total + qB_Total - QLoad_Total;

    % myfprintf(logging, fid, "Time Period = %d, macroItr = %d and Area = %d: " + ...
    %     "Projected savings in substation power flow by using batteries: " + ...
    %     "%f percent.\n", t, macroItr, Area, percentageSavings); % always true

    macroIterationPLosses(macroItr, Area) = PLoss;
    macroIterationQLosses(macroItr, Area) = macroIterationQLoss;
    % macroIterationQLosses(macroItr, Area) = QLoss;
    macroIterationPSaves(macroItr, Area) = percentageSavings;

     if fileOpenedFlag
        fclose(fid);
    end

end