function x_NoLoss = singlephaselin(busDataTable_pu_Area, branchDataTable_Area, v2_parent_Area, S_connection_Area, B0Vals_Area, isLeaf_Area, ...
    Area, numAreas, graphDFS_Area_Table, R_Area_Matrix, X_Area_Matrix, timePeriodNum, itr, varargin)

 % Default values for optional arguments
    verbose = false;
    logging = false;
    CVR = [0; 0];
    V_max = 1.05;
    V_min = 0.95;
    Qref_DER = 0.00;
    Vref_DER = 1.00;
    saveToFile = true;
    fileExtension = ".txt";
    systemName = "ieee123";
    saveLocationName = "logfiles/";
    delta_t = 0.25;
    etta_C = 0.80;
    etta_D = 0.80;
    chargeToPowerRatio = 4;
    soc_min = 0.30;
    soc_max = 0.95;
    % Process optional arguments
    numArgs = numel(varargin);

    if mod(numArgs, 2) ~= 0
        error('Optional arguments must be specified as name-value pairs.');
    end
    
    validArgs = ["verbose", "logging", "CVR", "V_max", "V_min", "saveToFile", "saveLocation", "Qref_DER", "Vref_DER", "fileExtension", "systemName", "delta_t", "etta_C", "etta_D", "chargeToPowerRatio", "soc_min", "soc_max"];
    
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
            case "fileExtension"
                fileExtension = argValue;
            case "systemName"
                fileExtension = argValue;
            case "delta_t"
                delta_t = argValue;
            case "etta_C"
                etta_C = argValue;
            case "etta_D"
                etta_D = argValue;
            case "chargeToPowerRatio"
                chargeToPowerRatio = argValue;
            case "soc_min"
                soc_min = argValue;
            case "soc_max"
                soc_max = argValue;
        end
    end
    
    strArea = convert2doubleDigits(Area);
    saveLocationFilename = strcat(saveLocationName, systemName, "/numAreas_", num2str(numAreas), "/Aeq_beq_area", strArea, "_singlephaselin", fileExtension);
    
    fileOpenedFlag = false;
    
    if itr ~= 1 || timePeriodNum ~= 1
        logging = false;
    end

    if logging && timePeriodNum == 1 && saveToFile && itr == 1 && Area == 2
        fileOpenedFlag = true;
        fid = fopen(saveLocationFilename, 'w');  % Open file for writing
    else
        logging = false;
        fid = 1;
    end
    
    N_Area = length(busDataTable_pu_Area.bus);
    m_Area = length(branchDataTable_Area.fb);
    fb_Area = branchDataTable_Area.fb;
    tb_Area = branchDataTable_Area.tb;
    P_L_Area = busDataTable_pu_Area.P_L;
    Q_L_Area = busDataTable_pu_Area.Q_L;
    Q_C_Area = busDataTable_pu_Area.Q_C;
    P_der_Area = busDataTable_pu_Area.P_der;
    S_der_Area = busDataTable_pu_Area.S_der;
    S_battMax_Area = S_der_Area;
    P_battMax_Area = P_der_Area;
    
    Emax_batt_Area = chargeToPowerRatio.*P_battMax_Area;

    busesWithDERs_Area = find(S_der_Area);
    nDER_Area = length(busesWithDERs_Area);
    busesWithBatts_Area = find(S_battMax_Area);
    nBatt_Area = length(busesWithBatts_Area);
    
    S_onlyDERbuses_Area = S_der_Area(busesWithDERs_Area);   %in PU
    
    P_onlyDERbuses_Area = P_der_Area(busesWithDERs_Area);   %in PU

    S_onlyBattBusesMax_Area = S_battMax_Area(busesWithBatts_Area);
    
    P_onlyBattBusesMax_Area = P_battMax_Area(busesWithBatts_Area);
    
    E_onlyBattBusesMax_Area = Emax_batt_Area(busesWithBatts_Area);

    lb_Pc_onlyBattBuses_Area = zeros*ones(nBatt_Area, 1);
    ub_Pc_onlyBattBuses_Area = P_onlyBattBusesMax_Area;
    
    lb_Pd_onlyBattBuses_Area = zeros*ones(nBatt_Area, 1);
    ub_Pd_onlyBattBuses_Area = P_onlyBattBusesMax_Area;
    
    lb_B_onlyBattBuses_Area = soc_min.*E_onlyBattBusesMax_Area;
    ub_B_onlyBattBuses_Area = soc_max.*E_onlyBattBusesMax_Area;

    lb_qD_onlyDERbuses_Area = -sqrt( S_onlyDERbuses_Area.^2 - P_onlyDERbuses_Area.^2 );
    ub_qD_onlyDERbuses_Area = sqrt( S_onlyDERbuses_Area.^2 - P_onlyDERbuses_Area.^2 );
    
    lb_qB_onlyBattBuses_Area = -sqrt( S_onlyBattBusesMax_Area.^2 - P_onlyBattBusesMax_Area.^2);
    ub_qB_onlyBattBuses_Area = sqrt( S_onlyBattBusesMax_Area.^2 - P_onlyBattBusesMax_Area.^2);
    
    if ~isLeaf_Area
        myfprintf(logging, fid, "Area %d is NOT a leaf area, does have child areas.\n", Area);
        for j = 1:size(S_connection_Area, 1)
            [P_L_Area(end-j+1), Q_L_Area(end-j+1)] = deal(real(S_connection_Area(end-j+1)), imag(S_connection_Area(end-j+1)));
        end
    else
        myfprintf(logging, fid, "Area %d does NOT have any child areas.\n", Area);
    end
    
    % numVarsNoLoss = [m_Area, m_Area, N_Area, nDER_Area];
    numVarsNoLoss = [m_Area, m_Area, N_Area, nDER_Area, nBatt_Area, nBatt_Area, nBatt_Area, nBatt_Area];
    ranges_noLoss = generateRangesFromValues(numVarsNoLoss);

    indices_P_noLoss = ranges_noLoss{1};
    indices_Q_noLoss = ranges_noLoss{2};
    indices_vFull_noLoss = ranges_noLoss{3};
    indices_v_noLoss = indices_vFull_noLoss(2:end);
    indices_qD_noLoss = ranges_noLoss{4};
    indices_B_noLoss = ranges_noLoss{5};
    indices_Pc_noLoss = ranges_noLoss{6};
    indices_Pd_noLoss = ranges_noLoss{7};
    indices_qB_noLoss = ranges_noLoss{8};

    Table_Area = [graphDFS_Area_Table.fbus graphDFS_Area_Table.tbus indices_P_noLoss' indices_Q_noLoss' indices_v_noLoss'];  % creating Table for variables P, Q ,l, V
    % Table_Area_Table = array2table(Table_Area, 'VariableNames', {'fbus', 'tbus', 'indices_P', 'indices_Q', 'indices_v'});
    
    myfprintf(logging, fid, "**********" + ...
        "Constructing Aeq and beq for Area %d.\n" + ...
        "***********\n", Area);

    % Initialization-
    
    CVR_P = CVR(1);                %% this will make the loads as constant power load
    CVR_Q = CVR(2);                %% this will make the loads as constant power load

    % numLinOptEquations = 3*m_Area + 1;
    numLinOptEquationsBFM = 3*m_Area + 1;
    numLinOptEquationsBFM_Batt = numLinOptEquationsBFM + nBatt_Area;
    numLinOptEquations = numLinOptEquationsBFM_Batt;
    % numOptVarsNoLoss = 2*m_Area + N_Area + nDER_Area;
    numOptVarsBFM_NoLoss = 2*m_Area + N_Area;
    numVarsBFM_DERs_NoLoss = numOptVarsBFM_NoLoss + nDER_Area;
    numOptVarsBFM_DERs_Batt_NoLoss = numVarsBFM_DERs_NoLoss + 4*nBatt_Area;
    numOptVarsNoLoss = numOptVarsBFM_DERs_Batt_NoLoss;
    Aeq_NoLoss = zeros(numLinOptEquations, numOptVarsNoLoss);
    beq_NoLoss = zeros(numLinOptEquations, 1);
    
    for currentBusNum = 2 : N_Area
        myfprintf(logging, fid, "*****\n" + ...
            "Checking for bus %d.\n" + ...
            "*****\n", currentBusNum);      

        parentBusIdx = find(tb_Area == currentBusNum) ;
        parentBusNum = fb_Area(parentBusIdx);
        
        myfprintf(logging, fid, "The parent of bus %d is bus %d at index %d.\n", currentBusNum, parentBusNum, parentBusIdx);

        PIdx = parentBusIdx;
        Aeq_NoLoss( PIdx, indices_P_noLoss(parentBusIdx) ) = 1;
        Aeq_NoLoss( PIdx, indices_v_noLoss(parentBusIdx) ) = -0.5 * CVR_P * P_L_Area( currentBusNum );
        
        QIdx = PIdx + m_Area;
        Aeq_NoLoss( QIdx, indices_Q_noLoss(parentBusIdx) ) = 1;
        Aeq_NoLoss( QIdx, indices_v_noLoss(parentBusIdx) ) = -0.5 * CVR_Q * Q_L_Area( currentBusNum );
        
        % List of Row Indices showing the set of 'children' buses 'under' our currentBus:
        childBusIndices = find(fb_Area == currentBusNum);
        if ~isempty(childBusIndices)
            Aeq_NoLoss(PIdx, indices_P_noLoss(childBusIndices) ) = -1;   % for P
            Aeq_NoLoss(QIdx, indices_Q_noLoss(childBusIndices) ) = -1;   % for Q
        end
        
        myfprintf(logging, fid, "Aeq(%d, P(%d)) = 1.\n", PIdx, parentBusIdx);
        for i = 1:length(childBusIndices)
            myfprintf(logging, fid, "Aeq(%d, P(%d)) = -1\n", PIdx, childBusIndices(i));
        end
        if CVR_P
            myfprintf(logging, fid, "Aeq(%d, v(%d)) = -0.5 * CVR_P * P_L(%d).\n", PIdx, parentBusIdx, currentBusNum);
        end
        
        myfprintf(logging, fid, "Aeq(%d, Q(%d)) = 1.\n", QIdx, parentBusIdx);
        for i = 1:length(childBusIndices)
            myfprintf(logging, fid, "Aeq(%d, Q(%d)) = -1\n", QIdx, childBusIndices(i));
        end
        if CVR_Q
            myfprintf(logging, fid, "Aeq(%d, v(%d)) = -0.5 * CVR_Q * Q_L(%d).\n", QIdx, parentBusIdx, currentBusNum);
        end

        % V equations
        vIdx = QIdx + m_Area;
        Aeq_NoLoss( vIdx, indices_v_noLoss(parentBusIdx) ) = 1;
        myfprintf(logging, fid, "Aeq(%d, v(%d)) = 1\n", vIdx, parentBusIdx);
    
        %Return the rows with the list of 'children' buses of 'under' the PARENT of our currentBus:
        %our currentBus will obviously also be included in the list.
        siblingBusesIndices = find(fb_Area == parentBusNum);
        siblingBuses = tb_Area(siblingBusesIndices);
        
        myfprintf(logging, fid, "The siblings of bus %d\n", currentBusNum);
        myfprintf(logging, fid, "include these buses: %d\n", siblingBuses)
        myfprintf(logging, fid, "at indices %d.\n", siblingBusesIndices);
        eldestSiblingIdx = siblingBusesIndices(1);
        eldestSiblingBus = siblingBuses(1);
        myfprintf(logging, fid,  "which makes bus %d at index %d as the eldest sibling.\n", eldestSiblingBus, eldestSiblingIdx);
        Aeq_NoLoss( vIdx, indices_vFull_noLoss( eldestSiblingIdx ) ) = -1;
        myfprintf(logging, fid, "Aeq(%d, v_Full(%d)) = -1\n", vIdx, eldestSiblingIdx);
        Aeq_NoLoss( vIdx, indices_P_noLoss(parentBusIdx) ) = 2 * R_Area_Matrix( parentBusNum, currentBusNum );
        myfprintf(logging, fid, "Aeq(%d, P(%d)) = 2*r(%d, %d).\n", vIdx, parentBusIdx, parentBusNum, currentBusNum);
        Aeq_NoLoss( vIdx, indices_Q_noLoss(parentBusIdx) ) = 2 * X_Area_Matrix( parentBusNum, currentBusNum );
        myfprintf(logging, fid, "Aeq(%d, Q(%d)) = 2*x(%d, %d).\n", vIdx, parentBusIdx, parentBusNum, currentBusNum);    
        
        beq_NoLoss(PIdx) = ...
            ( 1- 0.5 * CVR_P ) * ...
            ( P_L_Area( currentBusNum ) - P_der_Area( currentBusNum ) );
        myfprintf(logging, fid, "beq(%d) = (1 - 0.5*CVR_P)*(P_L(%d) - P_der(%d))\n", PIdx, currentBusNum, currentBusNum);
    
        beq_NoLoss(QIdx) =  ...
            ( 1- 0.5*CVR_Q ) * ...
            ( Q_L_Area( currentBusNum ) - Q_C_Area( currentBusNum ) );
        myfprintf(logging, fid, "beq(%d) = (1 - 0.5*CVR_Q)*(Q_L(%d) - Q_C(%d))\n", QIdx, currentBusNum, currentBusNum);
        
    end
    
    % substation voltage equation
    myfprintf(logging, fid, "And who can forget the substation voltage equation..\n")
    vSubIdx = 3*m_Area + 1;
    Aeq_NoLoss(vSubIdx, indices_vFull_noLoss(1) ) = 1;
    myfprintf(logging, fid, "Aeq(%d, v_Full(1)) = 1\n", vSubIdx);

    beq_NoLoss(vSubIdx) = v2_parent_Area;
    myfprintf(logging, fid, "beq(%d) = %.3f\n", vSubIdx, v2_parent_Area);

    Table_DER = zeros(nDER_Area, 5);
    
    myfprintf(logging, fid, "Now adding the required %d coefficients to model the %d DERs in just as many equations from among the %d branch reactive flow equations.\n", nDER_Area, nDER_Area, m_Area);
    for i = 1:nDER_Area
        currentBusNum = busesWithDERs_Area(i);
        parentBusIdx = find(graphDFS_Area_Table.tbus == currentBusNum);
        QIdx = parentBusIdx + m_Area;
        qD_Idx = indices_qD_noLoss(i);

        Aeq_NoLoss(QIdx, qD_Idx) = 1;
        myfprintf(logging, fid, "Aeq(%d, qD(%d)) = 1\n", QIdx, i);
        
        %setting other parameters for DGs:
        Table_DER(i, 2) = qD_Idx;
        
        % slope kq definiton:
        Table_DER(i, 3) = 2*ub_qD_onlyDERbuses_Area(i)/(V_max-V_min); % Qmax at Vmin, and vice versa
        
        % Q_ref, V_ref definition:
        Table_DER(i, 4) = Qref_DER;  %Qref
        Table_DER(i, 5) = Vref_DER;  %Vref
    end
    
    % Battery equation addition    
    myfprintf(logging, fid, "Now adding the required %d + %d = %d coefficients in just as many equations from among existing %d real and %d reactive branch flow equations, as well as %d new equations to model the SOCs for the batteries.\n", nBatt_Area, nBatt_Area, 2*nBatt_Area, m_Area, m_Area, nBatt_Area);
    for i = 1:nBatt_Area 
        currentBusNum = busesWithBatts_Area(i);
        myfprintf(logging, fid, "Battery number %d is actually bus number %d.\n", i, currentBusNum);
        parentBusIdx = find(tb_Area == currentBusNum);
        parentBusNum = fb_Area(parentBusIdx);
        myfprintf(logging, fid, "Whose parent is bus %d at index %d.\n", parentBusNum, parentBusIdx);
        PEqnIdx = parentBusIdx;
        QEqnIdx = parentBusIdx + m_Area;
        BEqnIdx = numLinOptEquationsBFM + i;
        myfprintf(logging, fid, "Obviously, this additional battery will go into equation #%d.\n", BEqnIdx);
        B_Idx = indices_B_noLoss(i);
        Pc_Idx = indices_Pc_noLoss(i);
        Pd_Idx = indices_Pd_noLoss(i);
        qB_Idx = indices_qB_noLoss(i);

        Aeq_NoLoss(PEqnIdx, Pc_Idx) = -1;
        myfprintf(logging, fid, "Aeq(%d, Pc(%d)) = -1\n", PEqnIdx, i);

        Aeq_NoLoss(PEqnIdx, Pd_Idx) = 1;
        myfprintf(logging, fid, "Aeq(%d, Pd(%d)) = 1\n", PEqnIdx, i);

        Aeq_NoLoss(QEqnIdx, qB_Idx) = 1;
        myfprintf(logging, fid, "Aeq(%d, qB(%d)) = 1\n", QEqnIdx, i);
        
        Aeq_NoLoss(BEqnIdx, B_Idx) = 1;
        myfprintf(logging, fid, "Aeq(%d, B(%d)) = 1\n", BEqnIdx, i);
        Aeq_NoLoss(BEqnIdx, Pc_Idx) = -delta_t*etta_C;
        myfprintf(logging, fid, "Aeq(%d, Pc(%d)) = -delta_t*etta_C\n", BEqnIdx, i);
        Aeq_NoLoss(BEqnIdx, Pd_Idx) = delta_t/etta_D;
        myfprintf(logging, fid, "Aeq(%d, Pd(%d)) = delta_t/etta_D\n", BEqnIdx, i);

        beq_NoLoss(BEqnIdx) = B0Vals_Area(i);
        myfprintf(logging, fid, "beq(%d) = B0(%d) = %f\n", BEqnIdx, i, B0Vals_Area(i));
    end

    numVarsForBoundsNoLoss = [1, numVarsNoLoss(1) - 1, numVarsNoLoss(2:3)]; % qD limits are specific to each machine, will be appended later.
    lbVals = [0, -5, -4, V_min^2];
    ubVals = [5, 5, 4, V_max^2];
    [lbBFM_NoLoss, ubBFM_NoLoss] = constructBoundVectors(numVarsForBoundsNoLoss, lbVals, ubVals);
    
    lbBFM_DER_NoLoss = [lbBFM_NoLoss; lb_qD_onlyDERbuses_Area];
    ubBFM_DER_NoLoss = [ubBFM_NoLoss; ub_qD_onlyDERbuses_Area];
    lbBFM_DER_Batt_NoLoss = [lbBFM_DER_NoLoss; lb_B_onlyBattBuses_Area; lb_Pc_onlyBattBuses_Area; lb_Pd_onlyBattBuses_Area; lb_qB_onlyBattBuses_Area];
    lb_NoLoss = lbBFM_DER_Batt_NoLoss;
    
    ubBFM_DER_Batt_NoLoss = [ubBFM_DER_NoLoss; ub_B_onlyBattBuses_Area; ub_Pc_onlyBattBuses_Area; ub_Pd_onlyBattBuses_Area; ub_qB_onlyBattBuses_Area];
    ub_NoLoss = ubBFM_DER_Batt_NoLoss;

    if fileOpenedFlag
        fclose(fid);
    end
        
    fBFM = zeros(numVarsBFM_DERs_NoLoss, 1);
    
    AeqBFM_NoLoss = Aeq_NoLoss(1:numLinOptEquationsBFM, 1:numVarsBFM_DERs_NoLoss);
    beqBFM_NoLoss = beq_NoLoss(1:numLinOptEquationsBFM);

    options = optimoptions('intlinprog','Display','off');
    
    myfprintf(verbose, "TimePeriod = %d, Area = %d and itr = %d: Performing Initialization Phase 1.\n", timePeriodNum, Area, itr);
    [xBFM_DER_NoLoss, ~, ~, ~] = intlinprog(fBFM, [], [], [], AeqBFM_NoLoss, beqBFM_NoLoss, lbBFM_DER_NoLoss, ubBFM_DER_NoLoss, options);
    
    if ~isempty(xBFM_DER_NoLoss)
        myfprintf(verbose, "TimePeriod = %d, Area = %d and itr = %d: Lossless Initialization WITHOUT batteries accomplished.\n", timePeriodNum, Area, itr);
    end
    
    P0_BFM_DER_NoLoss = xBFM_DER_NoLoss(indices_P_noLoss);
    Q0_BFM_DER_NoLoss = xBFM_DER_NoLoss(indices_Q_noLoss);
    v0_BFM_DER_NoLoss =  xBFM_DER_NoLoss(indices_vFull_noLoss);
    qD0_BFM_DER_NoLoss = xBFM_DER_NoLoss(indices_qD_noLoss);

    B0 = B0Vals_Area;
    Pc0 = zeros(nBatt_Area, 1);
    Pd0 = zeros(nBatt_Area, 1);
    qB0 = zeros(nBatt_Area, 1);
    
    x0_NoLoss = [P0_BFM_DER_NoLoss; Q0_BFM_DER_NoLoss; v0_BFM_DER_NoLoss; qD0_BFM_DER_NoLoss; B0; Pc0; Pd0; qB0];
    % [lb_NoLoss x0_NoLoss ub_NoLoss]
    flaggedForLimitViolation = false;
    for varNum = 1:numOptVarsNoLoss
        lbVal = lb_NoLoss(varNum);
        ubVal = ub_NoLoss(varNum);
        x0Val = x0_NoLoss(varNum);
        if lb_NoLoss(varNum) > x0_NoLoss(varNum)
            myfprintf(verbose, "TimePeriod = %d, Area = %d and itr = %d:  Oh no! x0_NoLoss(%d) < lb(%d) as %f < %f.\n", timePeriodNum, Area, itr, varNum, varNum, x0Val, lbVal);
            flaggedForLimitViolation = true;
        end
        if ub_NoLoss(varNum) < x0_NoLoss(varNum)
            myfprintf(verbose, "TimePeriod = %d, Area = %d and itr = %d:  Oh no! x0_NoLoss(%d) > ub(%d) as %f > %f.\n", timePeriodNum, Area, itr, varNum, varNum, x0Val, ubVal);
            flaggedForLimitViolation = true;
        end
    end
    
    if checkOptimalSolutionWithinBounds(x0_NoLoss, lb_NoLoss, ub_NoLoss)
        myfprintf("My native bound checker says that bounds are being violated.\n");
        error("Nani?");
    elseif flaggedForLimitViolation
        myfprintf(verbose, "TimePeriod = %d, Area = %d and itr = %d: x0_NoLoss within limits anyway? More like MATLAB stupid? Phase 1 of Initialization successful. Proceeding to Phase 2 of Initialization.\n", timePeriodNum, Area, itr)
    else
        myfprintf(verbose, "TimePeriod = %d, Area = %d and itr = %d: x0_NoLoss within limits. Phase 1 of Initialization successful. Proceeding to Phase 2 of Initialization.\n", timePeriodNum, Area, itr);
    end

    options = optimoptions('fmincon', 'Display', 'off', 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp');
    
    myfprintf(verbose, "TimePeriod = %d, Area = %d and itr = %d: Performing Initialization Phase 2\n", timePeriodNum, Area, itr);

    [x_NoLoss, ~, ~, ~] = ...
        fmincon(@(x)objfun(x, N_Area, nDER_Area, nBatt_Area, fb_Area, tb_Area, R_Area_Matrix, 'mainObjFun', "loss_min-fake", 'secondObjFun', "none", 'indices_Pd', indices_Pd_noLoss, 'indices_Pc', indices_Pc_noLoss), ...
        x0_NoLoss, [], [], Aeq_NoLoss, beq_NoLoss, lb_NoLoss, ub_NoLoss, [], options);
    
    if ~isempty(x_NoLoss)
        myfprintf(verbose, "TimePeriod = %d, Area = %d and itr = %d: Lossless Initialization WITH batteries accomplished.\n", timePeriodNum, Area, itr);
    else
        error("TimePeriod = %d, Area = %d and itr = %d: Lossless Initialization WITH batteries failed.\n", timePeriodNum, Area, itr);
    end
    
    % mydisplay(verbose, x_NoLoss)

end