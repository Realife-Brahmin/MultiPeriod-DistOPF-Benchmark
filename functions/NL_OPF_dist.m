function [v2_Area, S_Area, qD_Full_Area, BVals_Area,...
    microIterationLosses, itr, ...
    time_dist, R_Area_Matrix, graphDFS_Area, N_Area, m_Area, nBatt_Area, busDataTable_pu_Area, ...
    branchDataTable_Area] = ...
    ...
    NL_OPF_dist(v2_parent_Area, S_connection_Area, B0Vals_Area, ...
    Area, isLeaf_Area, isRoot_Area, numChildAreas_Area, numAreas, ...
    microIterationLosses, time_dist, itr, ...
    CB_FullTable, varargin)
    
 % Default values for optional arguments
    verbose = false;
    CVR = [0; 0];
    V_max = 1.05;
    V_min = 0.95;
    Qref_DER = 0.00;
    Vref_DER = 1.00;
    delta_t = 0.25;
    etta_D = 0.80;
    etta_C = 0.80;
    chargeToPowerRatio = 4;
    soc_min = 0.30;
    soc_max = 0.95;

    saveToFile = false;
    strArea = convert2doubleDigits(Area);
    fileExtension = ".txt";
    systemName = "ieee123";
    saveLocationName = "logfiles/";
    fileOpenedFlag = false;

    % Process optional arguments
    numArgs = numel(varargin);

    if mod(numArgs, 2) ~= 0
        error('Optional arguments must be specified as name-value pairs.');
    end
    
    validArgs = ["verbose", "systemName", "CVR", "V_max", "V_min", "saveToFile", "saveLocation", "Qref_DER", "Vref_DER", "fileExtension", "delta_t", "etta_C", "etta_D", "chargeToPowerRatio", "soc_min", "soc_max"];
    
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
    
    saveLocationFilename = strcat(saveLocationName , systemName, "/numAreas_", num2str(numAreas), "/Aeq_beq_area", strArea, fileExtension);
    if itr ~= 0 
        verbose = false;
    end

    if verbose && saveToFile && itr == 0 && Area == 2
        fileOpenedFlag = true;
        fid = fopen(saveLocationFilename, 'w');  % Open file for writing
    else
        verbose = false;
        fid = 1;
    end
    

    [busDataTable_pu_Area, branchDataTable_Area, edgeMatrix_Area, R_Area, X_Area] ...
        = extractAreaElectricalParameters(Area, itr, isRoot_Area, systemName, numAreas, CB_FullTable, numChildAreas_Area);
    
    N_Area = length(busDataTable_pu_Area.bus);
    m_Area = length(branchDataTable_Area.fb);
    fb_Area = branchDataTable_Area.fb;
    tb_Area = branchDataTable_Area.tb;
    P_L_Area = busDataTable_pu_Area.P_L
    Q_L_Area = busDataTable_pu_Area.Q_L;
    Q_C_Area = busDataTable_pu_Area.Q_C;
    P_der_Area = busDataTable_pu_Area.P_der
    S_der_Area = busDataTable_pu_Area.S_der;
    S_battMax_Area = S_der_Area;
    P_battMax_Area = P_der_Area;

    Emax_batt_Area = chargeToPowerRatio.*P_battMax_Area;

% Update the Parent Complex Power Vector with values from interconnection. 

    numChildAreas = size(S_connection_Area, 1); %could even be zero for a child-less area
    
    for j = 1:numChildAreas %nodes are numbered such that the last numChildAreas nodes are actually the interconnection nodes too.
        %The load of parent area at the node of interconnection is
        %basically the interconnection area power
        P_L_Area(end-j+1) = real(S_connection_Area(end-j+1));                 %in PU
        Q_L_Area(end-j+1) = imag(S_connection_Area(end-j+1));     
    end
    
    % DER Configuration:
    busesWithDERs_Area = find(S_der_Area); %all nnz element indices
    nDER_Area = length(busesWithDERs_Area);
    busesWithBatts_Area = find(S_battMax_Area);
    nBatt_Area = length(busesWithBatts_Area);
    B0Vals_Area = B0Vals_Area(1:nBatt_Area);

    myfprintf(verbose, fid, strcat("Number of DERs in Area ", num2str(Area), " : ", num2str(nDER_Area), ".\n") );
    myfprintf(verbose, fid, strcat("Number of Batteries in Area ", num2str(Area), " : ", num2str(nBatt_Area), ".\n") );

    S_onlyDERbuses_Area = S_der_Area(busesWithDERs_Area);   %in PU
    S_onlyBattBusesMax_Area = S_battMax_Area(busesWithBatts_Area);
    
    P_onlyDERbuses_Area = P_der_Area(busesWithDERs_Area);   %in PU
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

    graphDFS_Area = edgeMatrix_Area; %not doing any DFS
    graphDFS_Area_Table = array2table(graphDFS_Area, 'VariableNames', {'fbus', 'tbus'});

    
    R_Area_Matrix = zeros(N_Area, N_Area);
    X_Area_Matrix = zeros(N_Area, N_Area);
    
    % Matrix form of R and X in terms of graph
    for currentBusNum = 1: N_Area - 1
        R_Area_Matrix(fb_Area(currentBusNum), tb_Area(currentBusNum)) = R_Area(currentBusNum);
        R_Area_Matrix(tb_Area(currentBusNum), fb_Area(currentBusNum)) = R_Area_Matrix(fb_Area(currentBusNum), tb_Area(currentBusNum)) ;
        X_Area_Matrix(fb_Area(currentBusNum), tb_Area(currentBusNum)) = X_Area(currentBusNum);
        X_Area_Matrix(tb_Area(currentBusNum), fb_Area(currentBusNum)) = X_Area_Matrix(fb_Area(currentBusNum), tb_Area(currentBusNum)) ;
    end
% Initializing vectors to be used in the Optimization Problem formulation

    % defining the unknowns for phaseA
    
    % numVarsFull = [m_Area, m_Area, m_Area, N_Area, nDER_Area];
    numVarsFull = [m_Area, m_Area, m_Area, N_Area, nDER_Area, nBatt_Area, nBatt_Area, nBatt_Area, nBatt_Area];
    ranges_Full = generateRangesFromValues(numVarsFull);

    indices_P = ranges_Full{1};
    indices_Q = ranges_Full{2};
    indices_l = ranges_Full{3};
    indices_vFull = ranges_Full{4};
    indices_v = indices_vFull(2:end);
    indices_qD = ranges_Full{5};
    indices_B = ranges_Full{6};
    indices_Pc = ranges_Full{7};
    indices_Pd = ranges_Full{8};
    indices_qB = ranges_Full{9};

    % Table_Area = [fb_Area tb_Area indices_P' indices_Q' indices_l' indices_v'];  % creating Table for variables P, Q ,l, V
    % Table_Area = [fb_Area tb_Area indices_P' indices_Q' indices_l' indices_v' indices_qD' indices_B' indices_Pc' indices_Pd' indices_qB'];  % creating Table for variables P, Q ,l, V

    % Table_Area_Table = array2table(Table_Area, 'VariableNames', {'fbus', 'tbus', 'indices_P', 'indices_Q', 'indices_l', 'indices_v'});
    % Table_Area_Table = array2table(Table_Area, 'VariableNames', {'fbus', 'tbus', 'indices_P', 'indices_Q', 'indices_l', 'indices_v', 'indices_B', 'indices_Pc', 'indices_Pd', 'indices_qB'});

    % Initialization-
    
     myfprintf(verbose, fid, "**********" + ...
        "Constructing Aeq and beq for Area %d.\n" + ...
        "***********\n", Area);

    CVR_P = CVR(1);
    CVR_Q = CVR(2);
    
    % numLinOptEquations = 3*m_Area + 1;
    numLinOptEquations = 3*m_Area + 1 + nBatt_Area;
    % numOptVarsFull = 3*m_Area + N_Area + nDER_Area;
    numOptVarsFull = 3*m_Area + N_Area + nDER_Area + 4*nBatt_Area;
    Aeq = zeros(numLinOptEquations, numOptVarsFull);
    beq = zeros(numLinOptEquations, 1);

    for currentBusNum = 2 : N_Area
        myfprintf(verbose, fid, "*****\n" + ...
            "Checking for bus %d.\n" + ...
            "*****\n", currentBusNum);       
 
        % The row index showing the 'parent' bus of our currentBus:
        
        parentBusIdx = find(tb_Area == currentBusNum);
        parentBusNum = fb_Area(parentBusIdx);
        myfprintf(verbose, fid, "The parent of bus %d is bus %d at index %d.\n", currentBusNum, parentBusNum, parentBusIdx);

        PIdx = parentBusIdx;
        Aeq( PIdx, indices_P(parentBusIdx) ) = 1;
        Aeq( PIdx, indices_l(parentBusIdx) ) = -R_Area_Matrix( parentBusNum, currentBusNum );
        Aeq( PIdx, indices_v(parentBusIdx) ) = -0.5 * CVR_P * P_L_Area( currentBusNum );

        
        %Q equations
        QIdx = PIdx + m_Area;
        Aeq( QIdx, indices_Q(parentBusIdx) ) = 1;
        Aeq( QIdx, indices_l(parentBusIdx) ) = -X_Area_Matrix( parentBusNum, currentBusNum );
        Aeq( QIdx, indices_v(parentBusIdx) ) = -0.5 * CVR_Q * Q_L_Area( currentBusNum );

        
       % List of Row Indices showing the set of 'children' buses 'under' our currentBus:
        childBusIndices = find(fb_Area == currentBusNum);
        if ~isempty(childBusIndices)
            Aeq(PIdx, indices_P(childBusIndices) ) = -1;   % for P
            Aeq(QIdx, indices_Q(childBusIndices) ) = -1;   % for Q
        end
        
        myfprintf(verbose, fid, "Aeq(%d, P(%d)) = 1.\n", PIdx, parentBusIdx);
        myfprintf(verbose, fid, "Aeq(%d, l(%d)) = -r(%d, %d).\n", PIdx, parentBusIdx, parentBusNum, currentBusNum);
        for i = 1:length(childBusIndices)
            myfprintf(verbose, fid, "Aeq(%d, P(%d)) = -1\n", PIdx, childBusIndices(i));
        end
        if CVR_P
            myfprintf(verbose, fid, "Aeq(%d, v(%d)) = -0.5 * CVR_P * P_L(%d).\n", PIdx, parentBusIdx, currentBusNum);
        end
        

        myfprintf(verbose, fid, "Aeq(%d, Q(%d)) = 1.\n", QIdx, parentBusIdx);
        myfprintf(verbose, fid, "Aeq(%d, l(%d)) = -x(%d, %d).\n", QIdx, parentBusIdx, parentBusNum, currentBusNum);
        for i = 1:length(childBusIndices)
            myfprintf(verbose, fid, "Aeq(%d, Q(%d)) = -1\n", QIdx, childBusIndices(i));
        end
        if CVR_Q
            myfprintf(verbose, fid, "Aeq(%d, v(%d)) = -0.5 * CVR_Q * Q_L(%d).\n", QIdx, parentBusIdx, currentBusNum);
        end

        % V equations
        vIdx = QIdx + m_Area;
        Aeq( vIdx, indices_v(parentBusIdx) ) = 1;
        myfprintf(verbose, fid, "Aeq(%d, v(%d)) = 1\n", vIdx, parentBusIdx);

        %Return the rows with the list of 'children' buses of 'under' the PARENT of our currentBus:
        %our currentBus will obviously also be included in the list.
        siblingBusesIndices = find(fb_Area == parentBusNum);
        siblingBuses = tb_Area(siblingBusesIndices);

        myfprintf(verbose, fid, "The siblings of bus %d\n", currentBusNum);
        myfprintf(verbose, fid, "include these buses: %d\n", siblingBuses)
        myfprintf(verbose, fid, "at indices %d.\n", siblingBusesIndices);
        eldestSiblingIdx = siblingBusesIndices(1);
        eldestSiblingBus = siblingBuses(1);
        myfprintf(verbose, fid,  "which makes bus %d at index %d as the eldest sibling.\n", eldestSiblingBus, eldestSiblingIdx);
        Aeq( vIdx, indices_vFull( eldestSiblingIdx ) ) = -1;
        myfprintf(verbose, fid, "Aeq(%d, v_Full(%d)) = -1\n", vIdx, eldestSiblingIdx);
        Aeq( vIdx, indices_P(parentBusIdx) ) = 2 * R_Area_Matrix( parentBusNum, currentBusNum );
        myfprintf(verbose, fid, "Aeq(%d, P(%d)) = 2*r(%d, %d).\n", vIdx, parentBusIdx, parentBusNum, currentBusNum);
        Aeq( vIdx, indices_Q(parentBusIdx) ) = 2 * X_Area_Matrix( parentBusNum, currentBusNum );
        myfprintf(verbose, fid, "Aeq(%d, Q(%d)) = 2*x(%d, %d).\n", vIdx, parentBusIdx, parentBusNum, currentBusNum);
        Aeq( vIdx, indices_l(parentBusIdx) ) = ...
            -R_Area_Matrix( parentBusNum, currentBusNum )^2 + ...
            -X_Area_Matrix( parentBusNum, currentBusNum )^2 ;
        myfprintf(verbose, fid, "Aeq(%d, l(%d)) = -r(%d, %d)^2 -x(%d, %d)^2.\n", vIdx, parentBusIdx, parentBusNum, currentBusNum, parentBusNum, currentBusNum);
        

        beq(PIdx) = ...
            ( 1- 0.5 * CVR_P ) * ...
            ( P_L_Area( currentBusNum ) - P_der_Area( currentBusNum ) );
        myfprintf(verbose, fid, "beq(%d) = (1 - 0.5*CVR_P)*(P_L(%d) - P_der(%d))\n", PIdx, currentBusNum, currentBusNum);
    
        beq(QIdx) =  ...
            ( 1- 0.5*CVR_Q ) * ...
            ( Q_L_Area( currentBusNum ) - Q_C_Area( currentBusNum ) );
        myfprintf(verbose, fid, "beq(%d) = (1 - 0.5*CVR_Q)*(Q_L(%d) - Q_C(%d))\n", QIdx, currentBusNum, currentBusNum);

    end
    
    % substation voltage equation
    vSubIdx = 3*m_Area + 1;
    Aeq( vSubIdx, indices_vFull(1) ) = 1;
    myfprintf(verbose, fid, "Aeq(%d, v_Full(1)) = 1\n", vSubIdx);

    beq(vSubIdx) = v2_parent_Area;
    myfprintf(verbose, fid, "beq(%d) = %.3f\n", vSubIdx, v2_parent_Area);
    
    % DER equation addition
    Table_DER = zeros(nDER_Area, 5);
    
    for i = 1:nDER_Area
        currentBusNum = busesWithDERs_Area(i);
        parentBusIdx = find(tb_Area == currentBusNum);
        QIdx = parentBusIdx + m_Area;
        qD_Idx = indices_qD(i);
        Aeq(QIdx, qD_Idx) = 1;
        myfprintf(verbose, fid, "Aeq(%d, qD(%d)) = 1\n", QIdx, i);
        
        %setting other parameters for DGs:
        Table_DER(i, 2) = qD_Idx;
        
        % slope kq definiton:
        Table_DER(i, 3) = 2*ub_qD_onlyDERbuses_Area(i)/(V_max-V_min); % Qmax at Vmin, and vice versa
        
        % Q_ref, V_ref definition:
        Table_DER(i, 4) = Qref_DER;  %Qref
        Table_DER(i, 5) = Vref_DER;  %Vref
    end
    
    % Battery equation addition    
    for i = 1:nBatt_Area
        currentBusNum = busesWithBatts_Area(i);
        parentBusIdx = find(tb_Area == currentBusNum);
        PEqnIdx = parentBusIdx;
        QEqnIdx = parentBusIdx + m_Area;
        BEqnIdx = QEqnIdx + m_Area + N_Area;
        
        B_Idx = indices_B(i);
        Pc_Idx = indices_Pc(i);
        Pd_Idx = indices_Pd(i);
        qB_Idx = indices_qB(i);

        Aeq(PEqnIdx, Pc_Idx) = -1;
        myfprintf(verbose, fid, "Aeq(%d, Pc(%d)) = -1\n", PEqnIdx, i);

        Aeq(PEqnIdx, Pd_Idx) = 1;
        myfprintf(verbose, fid, "Aeq(%d, Pd(%d)) = 1\n", PEqnIdx, i);

        Aeq(QEqnIdx, qB_Idx) = 1;
        myfprintf(verbose, fid, "Aeq(%d, qB(%d)) = 1\n", QEqnIdx, i);
        
        Aeq(BEqnIdx, B_Idx) = 1;
        myfprintf(verbose, fid, "Aeq(%d, B(%d)) = 1\n", BEqnIdx, i);
        Aeq(BEqnIdx, Pc_Idx) = -delta_t*etta_C;
        myfprintf(verbose, fid, "Aeq(%d, Pc(%d)) = -delta_t*etta_C\n", BEqnIdx, i);
        Aeq(BEqnIdx, Pd_Idx) = delta_t/etta_D;
        myfprintf(verbose, fid, "Aeq(%d, Pd(%d)) = delta_t*etta_D\n", BEqnIdx, i);

        beq(BEqnIdx) = B0Vals_Area(i);
    end

    if fileOpenedFlag
        fclose(fid);
    end
    
    % calling linear solution for intial point
    x_linear_Area = singlephaselin(busDataTable_pu_Area, branchDataTable_Area, v2_parent_Area, S_connection_Area, B0Vals_Area, isLeaf_Area, ...
        Area, numAreas, graphDFS_Area_Table, R_Area_Matrix, X_Area_Matrix, itr, 'verbose', true);


    numVarsNoLoss = [m_Area, m_Area, N_Area, nDER_Area, nBatt_Area, nBatt_Area, nBatt_Area, nBatt_Area];
    ranges_noLoss = generateRangesFromValues(numVarsNoLoss);

    indices_P_noLoss = ranges_noLoss{1};
    indices_Q_noLoss = ranges_noLoss{2};
    indices_vFull_noLoss = ranges_noLoss{3};
    indices_qD_noLoss = ranges_noLoss{4};
    indices_B_noLoss = ranges_noLoss{5};
    indices_Pc_noLoss = ranges_noLoss{6};
    indices_Pd_noLoss = ranges_noLoss{7};
    indices_qB_noLoss = ranges_noLoss{8};
    
    P0_Area = x_linear_Area( indices_P_noLoss );
    Q0_Area = x_linear_Area( indices_Q_noLoss );
    v0_Area =  x_linear_Area( indices_vFull_noLoss );
    qD0_Area = x_linear_Area( indices_qD_noLoss );
    B0_Area = x_linear_Area(indices_B_noLoss);
    Pc0_Area = x_linear_Area(indices_Pc_noLoss);
    Pd0_Area = x_linear_Area(indices_Pd_noLoss);
    qB0_Area = x_linear_Area(indices_qB_noLoss);
    
    Iflow0_Area = zeros(m_Area, 1);

    for currentBusNum = 2 : N_Area
        parentBusIdx = find(tb_Area == currentBusNum);
        siblingBusesIndices = find(parentBusNum == fb_Area);
        Iflow0_Area( parentBusIdx ) = ( P0_Area(parentBusIdx)^2 + Q0_Area(parentBusIdx)^2 ) / v0_Area(siblingBusesIndices(1));
    end
    
    % x0_Area = [P0_Area; Q0_Area; Iflow0_Area; v0_Area; qD0_Area];
    x0_Area = [P0_Area; Q0_Area; Iflow0_Area; v0_Area; qD0_Area; B0_Area; Pc0_Area; Pd0_Area; qB0_Area];

    
    % Definig Limits
    
    numVarsForBoundsFull = [1, numVarsFull(1) - 1, numVarsFull(2:3) ]; % qD limits are specific to each machine, will be appended later.
    lbVals = [0, -1500, -1500, 0, V_min^2];
    ubVals = [1500, 1500, 1500, 1500, V_max^2];
    [lb_Area, ub_Area] = constructBoundVectors(numVarsForBoundsFull, lbVals, ubVals);
    
    % lb_AreaFull = [lb_Area; lb_qD_onlyDERbuses_Area];
    lb_AreaFull = [lb_Area; lb_qD_onlyDERbuses_Area; lb_B_onlyBattBuses_Area; lb_Pc_onlyBattBuses_Area; lb_Pd_onlyBattBuses_Area; lb_qB_onlyBattBuses_Area];

    % ub_AreaFull = [ub_Area; ub_qD_onlyDERbuses_Area];
    ub_AreaFull = [ub_Area; ub_qD_onlyDERbuses_Area; ub_B_onlyBattBuses_Area; ub_Pc_onlyBattBuses_Area; ub_Pd_onlyBattBuses_Area; ub_qB_onlyBattBuses_Area];

    
    if itr == 0 && Area == 2
        mydisplay(verbose, "branchTable",  graphDFS_Area_Table)
        mydisplay(verbose, "Aeq", Aeq)
        mydisplay(verbose, "beq", beq)
        mydisplay(verbose, "lb", lb_AreaFull)
        mydisplay(verbose, "ub", ub_AreaFull)
    end

    %  Optimization - 
    options = optimoptions('fmincon', 'Display', 'off', 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp');
    
    startSolvingForOptimization = tic;

    [x, ~, ~, ~] = fmincon( @(x)objfunTables(x, N_Area, fb_Area, tb_Area, indices_l, R_Area_Matrix), ...
                              x0_Area, [], [], Aeq, beq, lb_AreaFull, ub_AreaFull, ...
                              @(x)eqcons(x, Area, N_Area, ...
                              fb_Area, tb_Area, indices_P, indices_Q, indices_l, indices_vFull, ...
                              itr, systemName, numAreas, "verbose", false, "saveToFile", false),...
                              options);
    
    optimizationSolutionTime = toc(startSolvingForOptimization);
    
    time_dist(itr+1, Area) = optimizationSolutionTime;
    
    % Result
    P_Area = x(indices_P); %m_Areax1
    Q_Area = x(indices_Q); %m_Areax1
    S_Area = complex(P_Area, Q_Area); %m_Areax1
    v2_Area = x(indices_vFull); %N_Areax1
    v2_Area(1) = v2_parent_Area;
    
    qD_Area = x(indices_qD);

    qD_Full_Area = zeros(N_Area, 1);

    for i = 1 : nDER_Area
        qD_Full_Area( busesWithDERs_Area(i) ) = qD_Area(i);
    end
    
    BVals_Area = x(indices_B);
    
    B_Full_Area = zeros(N_Area, 1);

    for i = 1 : nBatt_Area
        B_Full_Area( busesWithBatts_Area(i) ) = BVals_Area(i);
    end

    microIterationLosses(itr + 1, Area) = P_Area(1) + sum(P_der_Area) - sum(P_L_Area);

end