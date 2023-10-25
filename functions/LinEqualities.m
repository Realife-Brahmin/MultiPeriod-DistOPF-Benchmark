function [Aeq, beq, lb_Area9, ub_Area9, x0, areaInfo] = LinEqualities(areaInfo, simInfo, lambdaVals, pvCoeffVals, v_parent_Area_1toT, varargin)

    % Default values
    % logging_Aeq_beq = false;
    T = simInfo.T;
    % Input parser
    p = inputParser;
    addParameter(p, 'delta_t', 0.25, @isnumeric);
    addParameter(p, 'etta_C', 0.95, @isnumeric);
    addParameter(p, 'etta_D', 0.95, @isnumeric);
    addParameter(p, 'V_max', 1.05, @isnumeric);
    addParameter(p, 'V_min', 0.95, @isnumeric);
    addParameter(p, 'verbose', false, @islogical);

    parse(p, varargin{:});
    
    kwargs = p.Results;
    delta_t = kwargs.delta_t;
    etta_C = kwargs.etta_C;
    etta_D = kwargs.etta_D;
    V_max = kwargs.V_max;
    V_min = kwargs.V_min;
    verbose = kwargs.verbose;

    % Unpacking areaInfo
    N_Area = areaInfo.N_Area;
    m_Area = areaInfo.m_Area;
    fb_Area = areaInfo.fb_Area;
    tb_Area = areaInfo.tb_Area;
    P_L_Area = areaInfo.P_L_Area;
    Q_L_Area = areaInfo.Q_L_Area;
    Q_C_Area = areaInfo.Q_C_Area;
    P_der_Area = areaInfo.P_der_Area;
    % S_der_Area = areaInfo.S_der_Area;
    % S_battMax_Area = areaInfo.S_battMax_Area;
    % P_battMax_Area = areaInfo.P_battMax_Area;
    % Emax_batt_Area = areaInfo.Emax_batt_Area;
    busesWithDERs_Area = areaInfo.busesWithDERs_Area;
    nDER_Area = areaInfo.nDER_Area;
    % busesWithBatts_Area = areaInfo.busesWithBatts_Area;
    nBatt_Area = areaInfo.nBatt_Area;
    % S_onlyDERbuses_Area = areaInfo.S_onlyDERbuses_Area;
    % P_onlyDERbuses_Area = areaInfo.P_onlyDERbuses_Area;
    % S_onlyBattBusesMax_Area = areaInfo.S_onlyBattBusesMax_Area;
    % P_onlyBattBusesMax_Area = areaInfo.P_onlyBattBusesMax_Area;
    % E_onlyBattBusesMax_Area = areaInfo.E_onlyBattBusesMax_Area;
    lb_Pc_onlyBattBuses_Area = areaInfo.lb_Pc_onlyBattBuses_Area;
    ub_Pc_onlyBattBuses_Area = areaInfo.ub_Pc_onlyBattBuses_Area;
    lb_Pd_onlyBattBuses_Area = areaInfo.lb_Pd_onlyBattBuses_Area;
    ub_Pd_onlyBattBuses_Area = areaInfo.ub_Pd_onlyBattBuses_Area;
    lb_B_onlyBattBuses_Area = areaInfo.lb_B_onlyBattBuses_Area;
    ub_B_onlyBattBuses_Area = areaInfo.ub_B_onlyBattBuses_Area;
    lb_qD_onlyDERbuses_Area = areaInfo.lb_qD_onlyDERbuses_Area;
    ub_qD_onlyDERbuses_Area = areaInfo.ub_qD_onlyDERbuses_Area;
    lb_qB_onlyBattBuses_Area = areaInfo.lb_qB_onlyBattBuses_Area;
    ub_qB_onlyBattBuses_Area = areaInfo.ub_qB_onlyBattBuses_Area;
    B0Vals_pu_Area = areaInfo.B0Vals_pu_Area;
    R_Area_Matrix = areaInfo.R_Area_Matrix;
    X_Area_Matrix = areaInfo.X_Area_Matrix;

    listNumVars1 = [m_Area*ones(1, 3), N_Area, nDER_Area, nBatt_Area*ones(1, 4)];
    nVars1 = sum(listNumVars1);
    nVarsT = nVars1 * T;
    listNumLinEqns1 = [m_Area*ones(1, 2), N_Area, nBatt_Area];
    nLinEqns1 = sum(listNumLinEqns1);
    nLinEqnsT = nLinEqns1 * T;
    
    Aeq = zeros(nLinEqnsT, nVarsT);
    sz = size(Aeq);
    beq = zeros(nLinEqnsT, 1);

    eqnIndicesT = generateRangesFromValuesT(listNumLinEqns1, T);
    varIndicesT = generateRangesFromValuesT(listNumVars1, T);
    % keyboard;
    

    indices_Pflow = eqnIndicesT{1};
    indices_Qflow = eqnIndicesT{2};
    indices_KVL = eqnIndicesT{3};
    indices_SOC = eqnIndicesT{4};

    indices_Pij = varIndicesT{1};
    indices_Qij = varIndicesT{2};
    indices_lij = varIndicesT{3};
    indices_vAllj = varIndicesT{4};
    indices_vj = excludeFirstElement(indices_vAllj);
    indices_qDj = varIndicesT{5};
    indices_Bj = varIndicesT{6};
    indices_Pdj = varIndicesT{7};
    indices_Pcj = varIndicesT{8};
    indices_qBj = varIndicesT{9};
    
    % Assuming you have already created the 'areaInfo' structure

    % Assign each variable as a field in areaInfo
    areaInfo.indices_Pflow = indices_Pflow;
    areaInfo.indices_Qflow = indices_Qflow;
    areaInfo.indices_KVL = indices_KVL;
    areaInfo.indices_SOC = indices_SOC;
    
    areaInfo.indices_Pij = indices_Pij;
    areaInfo.indices_Qij = indices_Qij;
    areaInfo.indices_lij = indices_lij;
    areaInfo.indices_vAllj = indices_vAllj;
    areaInfo.indices_vj = excludeFirstElement(indices_vAllj);  % Assuming excludeFirstElement is a defined function
    areaInfo.indices_qDj = indices_qDj;
    areaInfo.indices_Bj = indices_Bj;
    areaInfo.indices_Pdj = indices_Pdj;
    areaInfo.indices_Pcj = indices_Pcj;
    areaInfo.indices_qBj = indices_qBj;

    % [areaInfo.fb_Area areaInfo.tb_Area]
     for j = 2 : N_Area
        myfprintf(verbose, 1, "*****\n" + ...
            "Checking for bus %d.\n" + ...
            "*****\n", j);       
 
        % The row index showing the 'parent' bus of our currentBus:
        
        i_Idx = find(tb_Area == j);
        i = fb_Area(i_Idx);
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "The parent of bus %d is bus %d at index %d.\n", j, i, i_Idx);
        
        k_indices = find(fb_Area == j);

        
        js_indices = find(fb_Area==i);

        jes_Idx = js_indices(1);

        indices_Pflow_ij_T = getIndicesT(indices_Pflow, i_Idx);
        indices_Qflow_ij_T = getIndicesT(indices_Qflow, i_Idx);
        indices_KVL_ij_T = getIndicesT(indices_KVL, i_Idx);
        
        indices_Pij_T = getIndicesT(indices_Pij, i_Idx);
        indices_Qij_T = getIndicesT(indices_Qij, i_Idx);
        indices_lij_T = getIndicesT(indices_lij, i_Idx);
        indices_vj_T = getIndicesT(indices_vj, i_Idx);
        
        % PFlow equations in Aeq, beq | BFM variables only
        row = indices_Pflow_ij_T;
        Aeq(sub2ind(sz, row, indices_Pij_T)) = 1;
        for k_num = 1:length(k_indices)
            k_Idx = k_indices(k_num);
            indices_Pjk_T = getIndicesT(indices_Pij, k_Idx);
            Aeq(sub2ind(sz, row, indices_Pjk_T)) = -1;
        end
        Aeq(sub2ind(sz, row, indices_lij_T)) = -R_Area_Matrix(i, j);

        beq(row) = lambdaVals*P_L_Area(j);

        % QFlow equations in Aeq, beq | BFM variables only
        row = indices_Qflow_ij_T;

        Aeq(sub2ind(sz, row, indices_Qij_T)) = 1;
        for k_num = 1:length(k_indices)
            k_Idx = k_indices(k_num);
            indices_Qjk_T = getIndicesT(indices_Qij, k_Idx);
            Aeq(sub2ind(sz, row, indices_Qjk_T)) = -1;
        end
        Aeq(sub2ind(sz, row, indices_lij_T)) = -X_Area_Matrix(i, j);
        beq(row) = lambdaVals*Q_L_Area(j) - Q_C_Area(j);

        % Aeq_Full( QIdx, indices_v(i_Idx) ) = -0.5 * CVR_Q * Q_L_Area( j );

        % if CVR_P
        %     myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, v(%d)) = -0.5 * CVR_P * P_L(%d).\n", PIdx, i_Idx, j);
        % end
        
        % if CVR_Q
        %     myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, v(%d)) = -0.5 * CVR_Q * Q_L(%d).\n", QIdx, i_Idx, j);
        % end
        
        % KVL equations in Aeq, beq | Full equation (contains only BFM
        % variables)
        row = indices_KVL_ij_T;
        % Aeq(sub2ind(sz, row, indices_vj_T)) = 1;
        Aeq(sub2ind(sz, row, indices_vj_T)) = -1;

        indices_vi_T = getIndicesT(indices_vAllj, jes_Idx);
        % Aeq(sub2ind(sz, row, indices_vi_T)) = -1;
        Aeq(sub2ind(sz, row, indices_vi_T)) = 1;

        Aeq(sub2ind(sz, row, indices_lij_T)) = R_Area_Matrix(i, j)^2 + X_Area_Matrix(i, j)^2;
        Aeq(sub2ind(sz, row, indices_Pij_T)) = -2*R_Area_Matrix(i, j);
        Aeq(sub2ind(sz, row, indices_Qij_T)) = -2*X_Area_Matrix(i, j);

    end
    
    % 'Substation' Bus KVL | BFM variables required only
    indices_KVL_12_T = getIndicesT(indices_KVL, N_Area);
    indices_vAll1_T = getIndicesT(indices_vAllj, 1);
    row = indices_KVL_12_T;
    Aeq(sub2ind(sz, row, indices_vAll1_T)) = 1;
    beq(row) = v_parent_Area_1toT;
    
    for der_num = 1:nDER_Area
        j = busesWithDERs_Area(der_num);
        indices_qDj_T = getIndicesT(indices_qDj, der_num);
        
        row = indices_Pflow_ij_T;
        beq(row) = beq(row) - transpose(P_der_Area(j).*pvCoeffVals);

        row = indices_Qflow_ij_T;
        Aeq(sub2ind(sz, row, indices_qDj_T)) = 1;
        
    end
    
    for batt_num = 1:nBatt_Area
        
        indices_SOC_j_T = getIndicesT(indices_SOC, batt_num);
        indices_SOC_j_T_2toT = indices_SOC_j_T(2:T);
        indices_SOC_j_T_1 = indices_SOC_j_T(1);
        
        indices_Bj_T = getIndicesT(indices_Bj, batt_num);
        indices_Pdj_T = getIndicesT(indices_Pdj, batt_num);
        indices_Pcj_T = getIndicesT(indices_Pcj, batt_num);
        indices_qBj_T = getIndicesT(indices_qBj, batt_num);
        
        % Pflow equations | Battery Variables
        row = indices_Pflow_ij_T;
        Aeq(sub2ind(sz, row, indices_Pdj_T)) = 1;
        Aeq(sub2ind(sz, row, indices_Pcj_T)) = -1;
        
        % Qflow equations | Battery Variables
        row = indices_Qflow_ij_T;
        Aeq(sub2ind(sz, row, indices_qBj_T)) = 1;
        
        % SOC equations | Battery Variables | First Time Interval
        row = indices_SOC_j_T_1;
        Aeq(row, indices_Bj_T(1)) = -1;
        Aeq(row, indices_Pdj_T(1)) = delta_t*etta_C;
        Aeq(row, indices_Pcj_T(1)) = -delta_t/etta_D;
        beq(row) = -B0Vals_pu_Area(batt_num);

        % SOC equations | Battery Variables | First Time Interval
        row = indices_SOC_j_T_2toT;
        Aeq(sub2ind(sz, row, indices_Bj_T(2:T))) = -1;
        Aeq(sub2ind(sz, row, indices_Bj_T(1:T-1))) = 1;
        Aeq(sub2ind(sz, row, indices_Pcj_T(2:T))) = delta_t*etta_C;
        Aeq(sub2ind(sz, row, indices_Pdj_T(2:T))) = -delta_t/etta_D;

    end


    numVarsBFM4 = [1, listNumVars1(1) - 1, listNumVars1(2:4) ]; % qD limits are specific to each machine, will be appended later.
    lbVals4 = [0, -5, -15, 0, V_min^2];
    ubVals4 = [5, 5, 5, 15, V_max^2];
    [lb_Area4, ub_Area4] = constructBoundVectors(numVarsBFM4, lbVals4, ubVals4);
    lb_Area9 = repmat([lb_Area4; lb_qD_onlyDERbuses_Area; lb_B_onlyBattBuses_Area; lb_Pc_onlyBattBuses_Area; lb_Pd_onlyBattBuses_Area; lb_qB_onlyBattBuses_Area], T, 1);
    ub_Area9 = repmat([ub_Area4; ub_qD_onlyDERbuses_Area; ub_B_onlyBattBuses_Area; ub_Pc_onlyBattBuses_Area; ub_Pd_onlyBattBuses_Area; ub_qB_onlyBattBuses_Area], T, 1);

    x0 = repmat([zeros(3*m_Area, 1); 1.00*ones(N_Area, 1); zeros(nDER_Area, 1); B0Vals_pu_Area; zeros(3*nBatt_Area, 1)], T, 1);

end