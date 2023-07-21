function x_Area_Linear = singlephaselin(busDataTable_pu_Area, branchDataTable_Area, v2_parent_Area, S_connection_Area, isLeaf_Area, ...
    Area, numAreas, graphDFS_Area, graphDFS_Area_Table, R_Area_Matrix, X_Area_Matrix, ...
    lb_Q_onlyDERbuses_Area, ub_Q_onlyDERbuses_Area, itr, varargin)

 % Default values for optional arguments
    verbose = false;
    CVR = [0; 0];
    V_max = 1.05;
    V_min = 0.95;
    Qref_DER = 0.00;
    Vref_DER = 1.00;
    saveToFile = true;
    fileExtension = ".txt";
    systemName = "ieee123";
    saveLocationName = "logfiles/";
    
    % Process optional arguments
    numArgs = numel(varargin);

    if mod(numArgs, 2) ~= 0
        error('Optional arguments must be specified as name-value pairs.');
    end
    
    validArgs = ["verbose", "CVR", "V_max", "V_min", "saveToFile", "saveLocation", "Qref_DER", "Vref_DER", "fileExtension", "systemName"];
    
    for i = 1:2:numArgs
        argName = varargin{i};
        argValue = varargin{i+1};
        
        if ~ischar(argName) || ~any(argName == validArgs)
            error('Invalid optional argument name.');
        end
        
        switch argName
            case "verbose"
                verbose = argValue;
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
        end
    end
    
    strArea = convert2doubleDigits(Area);
    saveLocationFilename = strcat(saveLocationName, systemName, "/numAreas_", num2str(numAreas), "/Aeq_beq_area", strArea, "_singlephaselin", fileExtension);
    
    fileOpenedFlag = false;
    
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
    
    N_Area = length(busDataTable_pu_Area.bus);
    m_Area = length(branchDataTable_Area.fb);
    fb_Area = branchDataTable_Area.fb;
    tb_Area = branchDataTable_Area.tb;
    P_L_Area = busDataTable_pu_Area.P_L;
    Q_L_Area = busDataTable_pu_Area.Q_L;
    Q_C_Area = busDataTable_pu_Area.Q_C;
    P_der_Area = busDataTable_pu_Area.P_der;
    S_der_Area = busDataTable_pu_Area.S_der;
    busesWithDERs_Area = find(S_der_Area);
    nDER_Area = length(busesWithDERs_Area);

    if ~isLeaf_Area
        myfprintf(verbose, fid, "Area %d is NOT a leaf area, does have child areas.\n", Area);
        for j = 1:size(S_connection_Area, 1)
            [P_L_Area(end-j+1), Q_L_Area(end-j+1)] = deal(real(S_connection_Area(end-j+1)), imag(S_connection_Area(end-j+1)));
        end
    else
        myfprintf(verbose, fid, "Area %d does NOT have any child areas.\n", Area);
    end
    
    numVarsNoLoss = [m_Area, m_Area, N_Area, nDER_Area];

    ranges_noLoss = generateRangesFromValues(numVarsNoLoss);

    indices_P = ranges_noLoss{1};
    indices_Q = ranges_noLoss{2};
    indices_vFull = ranges_noLoss{3};
    indices_v = indices_vFull(2:end);
    indices_qD = ranges_noLoss{4};
    

    numVarsForBoundsNoLoss = [1, numVarsNoLoss(1) - 1, numVarsNoLoss(2:end-1) ]; % qD limits are specific to each machine, will be appended later.
    lbVals = [0, -1500, -1500, V_min^2];
    ubVals = [1500, 1500, 1500, V_max^2];
    [lb_Area, ub_Area] = constructBoundVectors(numVarsForBoundsNoLoss, lbVals, ubVals);

    Table_Area = [graphDFS_Area_Table.fbus graphDFS_Area_Table.tbus indices_P' indices_Q' indices_v'];  % creating Table for variables P, Q ,l, V
    Table_Area_Table = array2table(Table_Area, 'VariableNames', {'fbus', 'tbus', 'indices_P', 'indices_Q', 'indices_v'});
    
    myfprintf(verbose, fid, "**********" + ...
        "Constructing Aeq and beq for Area %d.\n" + ...
        "***********\n", Area);

    % Initialization-
    
    CVR_P = CVR(1);                %% this will make the loads as constant power load
    CVR_Q = CVR(2);                %% this will make the loads as constant power load

    numLinOptEquations = 3*m_Area + 1;
    numOptVarsFull = 2*m_Area + N_Area + nDER_Area;
    Aeq = zeros(numLinOptEquations, numOptVarsFull);
    beq = zeros(numLinOptEquations, 1);
    
    for currentBusNum = 2 : N_Area
        myfprintf(verbose, fid, "*****\n" + ...
            "Checking for bus %d.\n" + ...
            "*****\n", currentBusNum);      

        parentBusIdx = find(tb_Area == currentBusNum) ;
        parentBusNum = fb_Area(parentBusIdx);
        
        myfprintf(verbose, fid, "The parent of bus %d is bus %d at index %d.\n", currentBusNum, parentBusNum, parentBusIdx);

        PIdx = parentBusIdx;
        Aeq( PIdx, indices_P(parentBusIdx) ) = 1;
        Aeq( PIdx, indices_v(parentBusIdx) ) = -0.5 * CVR_P * P_L_Area( currentBusNum );
        
        QIdx = PIdx + m_Area;
        Aeq( QIdx, indices_Q(parentBusIdx) ) = 1;
        Aeq( QIdx, indices_v(parentBusIdx) ) = -0.5 * CVR_Q * Q_L_Area( currentBusNum );
        
        % List of Row Indices showing the set of 'children' buses 'under' our currentBus:
        childBusIndices = find(fb_Area == currentBusNum);
        if ~isempty(childBusIndices)
            Aeq(PIdx, indices_P(childBusIndices) ) = -1;   % for P
            Aeq(QIdx, indices_Q(childBusIndices) ) = -1;   % for Q
        end
        
        myfprintf(verbose, fid, "Aeq(%d, P(%d)) = 1.\n", PIdx, parentBusIdx);
        for i = 1:length(childBusIndices)
            myfprintf(verbose, fid, "Aeq(%d, P(%d)) = -1\n", PIdx, childBusIndices(i));
        end
        if CVR_P
            myfprintf(verbose, fid, "Aeq(%d, v(%d)) = -0.5 * CVR_P * P_L(%d).\n", PIdx, parentBusIdx, currentBusNum);
        end
        
        myfprintf(verbose, fid, "Aeq(%d, Q(%d)) = 1.\n", QIdx, parentBusIdx);
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

    Table_DER = zeros(nDER_Area, 5);
    
    for i = 1:nDER_Area
        currentBusNum = busesWithDERs_Area(i);
        parentBusIdx = find(graphDFS_Area_Table.tbus == currentBusNum);
        QIdx = parentBusIdx + m_Area;
        qD_Idx = indices_qD(i);

        Aeq(QIdx, qD_Idx) = 1;
        myfprintf(verbose, fid, "Aeq(%d, qD(%d)) = 1\n", QIdx, i);
        
        %setting other parameters for DGs:
        Table_DER(i, 2) = qD_Idx;
        
        % slope kq definiton:
        Table_DER(i, 3) = 2*ub_Q_onlyDERbuses_Area(i)/(V_max-V_min); % Qmax at Vmin, and vice versa
        
        % Q_ref, V_ref definition:
        Table_DER(i, 4) = Qref_DER;  %Qref
        Table_DER(i, 5) = Vref_DER;  %Vref
    end
    
    
    if fileOpenedFlag
        fclose(fid);
    end
    Tnvar = size(Aeq,2);         % total number of variables
        
    f = zeros(Tnvar,1);
    f(Table_Area(1,3)) = 0;
    
    lb_AreaFull = [lb_Area ;lb_Q_onlyDERbuses_Area];
    ub_AreaFull = [ub_Area; ub_Q_onlyDERbuses_Area];
    
    options = optimoptions('intlinprog','Display','off');
    [x_Area_Linear, ~, ~, ~] = intlinprog(f, [], [], [], Aeq, beq, lb_AreaFull, ub_AreaFull, options);

end