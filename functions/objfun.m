function [f, f_1toT] = objfun(x, simInfo, sysInfo, areaInfo, T, varargin)
    
    noBatteries = simInfo.noBatteries;
    % Default values for optional arguments
    verbose = false;
    % etta_C = 0.95;
    etta_C = simInfo.etta_C;
    % etta_D = 0.95;
    etta_D = simInfo.etta_D;
    % alpha = 3e-5;
    % alpha = 3e-4;
    % alpha = 3e-3;
    % alpha = 1e-3;
    alpha = simInfo.alpha;
    % alpha = 10;
    % alpha = 1;
    % alpha = 9e-4;
    % gamma = 1e3;
    % gamma = 10;
    % gamma = 1e0;
    gamma = simInfo.gamma;
    costArray = simInfo.costArray;
    delta_t = simInfo.delta_t;
    % mainObjFun = "func_PLoss";
    % secondObjFun = "func_SCD";
    % Process optional arguments
    kVA_B = sysInfo.kVA_B;
    kV_B = sysInfo.kV_B;
    % Default objective functions
    if ~noBatteries
        % objectiveFuns = {"func_PLoss", "func_SCD", "func_netChangeInSOC"}; 
        objectiveFuns = {"func_PSubsCost", "func_SCD", "func_netChangeInSOC"};
    else
        % objectiveFuns = {"func_PLoss"};
        objectiveFuns = {"func_PSubsCost"};
    end

    numArgs = numel(varargin);

    if mod(numArgs, 2) ~= 0
        error('Optional arguments must be specified as name-value pairs.');
    end
    
    % validArgs = ["verbose", "etta_C", "etta_D", "alpha", "gamma", "objectiveFuns"];
    validArgs = ["verbose", "objectiveFuns"];

    
    for i = 1:2:numArgs
        argName = varargin{i};
        argValue = varargin{i+1};
        
        if ~ischar(argName) || ~any(argName == validArgs)
            error('Invalid optional argument name.');
        end
        
        switch argName
            case "verbose"
                verbose = argValue;
            case "objectiveFuns"
                objectiveFuns = argValue;
        end
    end
    
    
    % Unpacking areaInfo
    N_Area = areaInfo.N_Area;
    m_Area = areaInfo.m_Area;
    fb_Area = areaInfo.fb_Area;
    tb_Area = areaInfo.tb_Area;
    % P_L_Area = areaInfo.P_L_Area;
    % Q_L_Area = areaInfo.Q_L_Area;
    % Q_C_Area = areaInfo.Q_C_Area;
    % P_der_Area = areaInfo.P_der_Area;
    % S_der_Area = areaInfo.S_der_Area;
    % S_battMax_Area = areaInfo.S_battMax_Area;
    % P_battMax_Area = areaInfo.P_battMax_Area;
    % Emax_batt_Area = areaInfo.Emax_batt_Area;
    % busesWithDERs_Area = areaInfo.busesWithDERs_Area;
    nDER_Area = areaInfo.nDER_Area;
    % busesWithBatts_Area = areaInfo.busesWithBatts_Area;
    if ~noBatteries
        nBatt_Area = areaInfo.nBatt_Area;
        B0Vals_pu_Area = areaInfo.B0Vals_pu_Area;
    else
        nBatt_Area = 0;
        B0Vals_pu_Area = [];
    end
    % S_onlyDERbuses_Area = areaInfo.S_onlyDERbuses_Area;
    % P_onlyDERbuses_Area = areaInfo.P_onlyDERbuses_Area;
    % S_onlyBattBusesMax_Area = areaInfo.S_onlyBattBusesMax_Area;
    % P_onlyBattBusesMax_Area = areaInfo.P_onlyBattBusesMax_Area;
    % E_onlyBattBusesMax_Area = areaInfo.E_onlyBattBusesMax_Area;
    % lb_Pc_onlyBattBuses_Area = areaInfo.lb_Pc_onlyBattBuses_Area;
    % ub_Pc_onlyBattBuses_Area = areaInfo.ub_Pc_onlyBattBuses_Area;
    % lb_Pd_onlyBattBuses_Area = areaInfo.lb_Pd_onlyBattBuses_Area;
    % ub_Pd_onlyBattBuses_Area = areaInfo.ub_Pd_onlyBattBuses_Area;
    % lb_B_onlyBattBuses_Area = areaInfo.lb_B_onlyBattBuses_Area;
    % ub_B_onlyBattBuses_Area = areaInfo.ub_B_onlyBattBuses_Area;
    % lb_qD_onlyDERbuses_Area = areaInfo.lb_qD_onlyDERbuses_Area;
    % ub_qD_onlyDERbuses_Area = areaInfo.ub_qD_onlyDERbuses_Area;
    % lb_qB_onlyBattBuses_Area = areaInfo.lb_qB_onlyBattBuses_Area;
    % ub_qB_onlyBattBuses_Area = areaInfo.ub_qB_onlyBattBuses_Area;
    R_Area_Matrix = areaInfo.R_Area_Matrix;
    X_Area_Matrix = areaInfo.X_Area_Matrix;

    listNumVars1 = [m_Area*ones(1, 3), N_Area, nDER_Area, nBatt_Area*ones(1, 4)];
    % nVars1 = sum(listNumVars1);
    % nVarsT = nVars1 * T;
    % listNumEqns1 = [m_Area*ones(1, 2), N_Area, nBatt_Area];
    % nEqns1 = sum(listNumEqns1);
    % nEqnsT = nEqns1 * T;
    
    % Aeq = zeros(nEqnsT, nVarsT);
    % sz = size(Aeq);
    % beq = zeros(nEqnsT, 1);

    % eqnIndicesT = generateRangesFromValuesT(listNumEqns1, T);
    varIndicesT = generateRangesFromValuesT(listNumVars1, T);
    

    % indices_Pflow = eqnIndicesT{1};
    % indices_Qflow = eqnIndicesT{2};
    % indices_KVL = eqnIndicesT{3};
    % indices_SOC = eqnIndicesT{4};
    
    indices_Pij = varIndicesT{1};

    % indices_Qij = varIndicesT{2};
    indices_lij = varIndicesT{3};
    % indices_vAllj = varIndicesT{4};
    % indices_vj = excludeFirstElement(indices_vAllj);
    % indices_qDj = varIndicesT{5};
    if ~noBatteries
        indices_Bj = varIndicesT{6};
        indices_Pdj = varIndicesT{7};
        indices_Pcj = varIndicesT{8};
    end
    % indices_qBj = varIndicesT{9};

    f = 0;
    f_1toT = zeros(T, 1);

    for objfun_num = 1:length(objectiveFuns)
        objFun = objectiveFuns{objfun_num};
        switch objFun
            case "func_PSubs"
                indices_PSubs_1toT = getIndicesT(indices_Pij, 1);
                row = indices_PSubs_1toT;
                f = f + sum(x(row)); % will give a value of cents
            
            case "func_PSubs_1toT" % don't use it for optimization
                for t = 1 : T
                    indices_PSubs_1toT = getIndicesT(indices_Pij, 1);
                    indices_PSubs_t = indices_PSubs_1toT(t);
                    row = indices_PSubs_t;
                    f_1toT(t) = f_1toT(t) + sum(x(row));
                end

            case "func_PSubsCost"
                indices_PSubs_1toT = getIndicesT(indices_Pij, 1);
                row = indices_PSubs_1toT;
                f = f + sum(x(row) .* costArray)*kVA_B*delta_t; % will give a value of cents
            
            case "func_PSubsCost_1toT" % don't use it for optimization
                for t = 1 : T
                    indices_PSubs_1toT = getIndicesT(indices_Pij, 1);
                    indices_PSubs_t = indices_PSubs_1toT(t);
                    row = indices_PSubs_t;
                    f_1toT(t) = f_1toT(t) + sum(x(row) .* costArray(t)) * kVA_B * delta_t;
                end

            case "func_PLoss"
                % ... your code for func_PLoss
                for j = 2 : N_Area
                    i_Idx = find(tb_Area == j);
                    i = fb_Area(i_Idx);

                    indices_lij_T = getIndicesT(indices_lij, i_Idx);
                    row = indices_lij_T;

                    if ~isempty(i_Idx)
                        f = f + sum(x(row) * R_Area_Matrix(i, j));
                    end
                end

            case "func_PLoss_1toT"

                f_1toT = zeros(T, 1);
                for t = 1 : T
                    for j = 2 : N_Area
                        i_Idx = find(tb_Area == j);
                        i = fb_Area(i_Idx);

                        indices_lij_1toT = getIndicesT(indices_lij, i_Idx);
                        indices_lij_t = indices_lij_1toT(t, :);
                        row = indices_lij_t;

                        if ~isempty(i_Idx)
                            f_1toT(t) = f_1toT(t) + sum(x(row) * R_Area_Matrix(i, j));
                        end
                    end
                end


            case "func_QLoss"
                % ... your code for func_QLoss
                for j = 2:N_Area
                    i_Idx = find(tb_Area == j);
                    i = fb_Area(i_Idx);
                    indices_lij_T = getIndicesT(indices_lij, i_Idx);
                    row = indices_lij_T;
                    f = f + sum(x(row) * X_Area_Matrix(i, j));
                end

            case "func_SCD"
                % ... your code for func_SCD
                myfprintf(verbose, "Complementarity Needs to Be Respected.\n");

                for batt_num = 1:nBatt_Area
                    indices_Pcj_T = getIndicesT(indices_Pcj, batt_num);
                    indices_Pdj_T = getIndicesT(indices_Pdj, batt_num);
                    f = f + alpha* ( sum(x(indices_Pcj_T)*(1-etta_C)) + sum(x(indices_Pdj_T)*(1/etta_D - 1)) );
                end
            
            case "func_netChangeInSOC"
                myfprintf(verbose, "SOCs need to be brought back to original at the end of horizon.\n");

                for batt_num = 1:nBatt_Area
                    indices_Bj_T = getIndicesT(indices_Bj, batt_num);
                    idx_Bj_T_T = indices_Bj_T(T);
                    Bj_T_0 = B0Vals_pu_Area(batt_num);
                    f = f + gamma* ( x(idx_Bj_T_T) - Bj_T_0 )^2;
                end

            otherwise
                error("Unknown objective function name: %s", objFun)

        end
    end

end