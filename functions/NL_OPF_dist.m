function [v2_Area, S_parent_Area, S_child_Area, ...
    microIterationLosses, decisionVars_Area, itr, ...
    time_dist, R_Area_Matrix, graphDFS_Area, N_Area, busDataTable_Area, ...
    branchDataTable_Area] = ...
    ...
    NL_OPF_dist(v2_parent_Area, S_connection_Area, ...
    Area, isLeaf_Area, isRoot_Area, numChildAreas_Area, N, systemName, numAreas, ...
    microIterationLosses, time_dist, itr, ...
    CB_FullTable, varargin)
    
 % Default values for optional arguments
    verbose = false;
    CVR = [0; 0];
    V_max = 1.05;
    V_min = 0.95;

    saveToFile = false;
    strArea = convert2doubleDigits(Area);
    fileExtension = ".txt";
    % systemName = "ieee123";
    saveLocationFilename = strcat("logfiles/", systemName, "/numAreas_", num2str(numAreas), "/Aeq_beq_area", strArea, fileExtension);
    fileOpenedFlag = false;

    % Process optional arguments
    numArgs = numel(varargin);

    if mod(numArgs, 2) ~= 0
        error('Optional arguments must be specified as name-value pairs.');
    end
    
    validArgs = ["verbose", "CVR", "V_max", "V_min", "saveToFile", "saveLocation"];
    
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
                saveLocationFilename = argValue;
        end
    end

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
% ExtractAreaElectricalParameters extracts electrical data required for OPF from an area's csv files. Optionally it can plot the graphs for the areas and save them as pngs.

    [busDataTable_Area, N_Area, branchDataTable_Area, fb_Area, tb_Area, edgeMatrix_Area, R_Area, X_Area, P_L_Area, Q_L_Area, Q_C_Area, P_der_Area, S_der_Area] ...
        = extractAreaElectricalParameters(Area, itr, N, isRoot_Area, systemName, numAreas, CB_FullTable, numChildAreas_Area);
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
    numBusesWithDERs_Area = length(busesWithDERs_Area);

    mydisp(verbose, ['Number of DERs in Area ', num2str(Area), ' : ', num2str(numBusesWithDERs_Area)]);
    
    S_onlyDERbuses_Area = S_der_Area(busesWithDERs_Area);   %in PU
    P_onlyDERbuses_Area = P_der_Area(busesWithDERs_Area);   %in PU
    lb_Q_onlyDERbuses_Area = -sqrt( S_onlyDERbuses_Area.^2 - P_onlyDERbuses_Area.^2 );
    ub_Q_onlyDERbuses_Area = sqrt( S_onlyDERbuses_Area.^2 - P_onlyDERbuses_Area.^2 );
    
    % graphDFS_Area = dfsearch(graph_Area, 1, 'edgetonew');
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
    %Assume N_Area = 41 (Area 1)
    m_Area = N_Area - 1;
    indices_P = 1:N_Area-1;  %1:40                                % defining the variables for P
    indices_Q = indices_P + (N_Area - 1); %41:80                             % defining the variables for Q
    indices_l = indices_Q + m_Area; %81:120                      % defining the variables for l(I^2)
    indices_v = indices_l + 1 + m_Area; %122:161                     % defining the variables for V
    Table_Area = [graphDFS_Area_Table.fbus graphDFS_Area_Table.tbus indices_P' indices_Q' indices_l' indices_v'];  % creating Table for variables P, Q ,l, V
    Table_Area_Table = array2table(Table_Area, 'VariableNames', {'fbus', 'tbus', 'indices_P', 'indices_Q', 'indices_l', 'indices_v'});
    indices_v_Full = transpose( indices_v(1)-1:indices_v(end) ) ; %121:161                   % voltage variables including parent node.

    % Initialization-
    
     myfprintf(verbose, fid, "**********" + ...
        "Constructing Aeq and beq for Area %d.\n" + ...
        "***********\n", Area);

    CVR_P = CVR(1);
    CVR_Q = CVR(2);
    
    % A and b matrix formulation-
    % mydisplay(verbose, graphDFS_Area_Table)
    % mydisplay(verbose, Table_Area_Table)
    
    for currentBusNum = 2 : N_Area
        myfprintf(verbose, fid, "*****\n" + ...
            "Checking for bus %d.\n" + ...
            "*****\n", currentBusNum);       
 
        % The row index showing the 'parent' bus of our currentBus:
        
        parentBusIdx = find(graphDFS_Area_Table.tbus == currentBusNum);
        parentBusNum = graphDFS_Area_Table.fbus(parentBusIdx);
        myfprintf(verbose, fid, "The parent of bus %d is bus %d at index %d.\n", currentBusNum, parentBusNum, parentBusIdx);

        % Aeq = zeros( 3*(N_Area-1), Table_Area_Table{end, end} ); %zeros(120, 121)
        % beq = zeros( 3*(N_Area-1), 1); %zeros(120, 1)
        % Aeq formulations
        %P equations
    %     busesWithDERs_Area
        PIdx = parentBusIdx;
        Aeq( PIdx, indices_P(parentBusIdx) ) = 1;
        myfprintf(verbose, fid, "Aeq(%d, P(%d)) = 1.\n", PIdx, parentBusIdx);
        Aeq( PIdx, indices_l(parentBusIdx) ) = -R_Area_Matrix( parentBusNum, currentBusNum );
        myfprintf(verbose, fid, "Aeq(%d, l(%d)) = -r(%d, %d).\n", PIdx, parentBusIdx, parentBusNum, currentBusNum);
        Aeq( PIdx, indices_v(parentBusIdx) ) = -0.5 * CVR_P * P_L_Area( currentBusNum );
        if CVR_P
            myfprintf(verbose, fid, "Aeq(%d, v(%d)) = -0.5 * CVR_P * P_L(%d).\n", PIdx, parentBusIdx, currentBusNum);
        end
        
        %Q equations
        QIdx = PIdx + (N_Area-1);
        Aeq( QIdx, indices_Q(parentBusIdx) ) = 1;
        myfprintf(verbose, fid, "Aeq(%d, Q(%d)) = 1.\n", QIdx, parentBusIdx);
        myfprintf(verbose, fid, "Note that we've used QIdx = %d, instead of %d for indexing into Aeq.\n", QIdx, indices_Q(parentBusIdx));
        Aeq( QIdx, indices_l(parentBusIdx) ) = -X_Area_Matrix( parentBusNum, currentBusNum );
        myfprintf(verbose, fid, "Aeq(%d, l(%d)) = -x(%d, %d).\n", QIdx, parentBusIdx, parentBusNum, currentBusNum);
        Aeq( QIdx, indices_v(parentBusIdx) ) = -0.5 * CVR_Q * Q_L_Area( currentBusNum );
        if CVR_Q
            myfprintf(verbose, fid, "Aeq(%d, v(%d)) = -0.5 * CVR_Q * Q_L(%d).\n", QIdx, parentBusIdx, currentBusNum);
        end
        
       % List of Row Indices showing the set of 'children' buses 'under' our currentBus:
        childBusIndices = find(graphDFS_Area_Table.fbus == currentBusNum);
        childBuses = graphDFS_Area_Table.tbus(childBusIndices);
        if isempty(childBusIndices)
            myfprintf(verbose, fid, "It is a leaf node.\n");
        else
            myfprintf(verbose, fid, "The child buses of bus %d\n", currentBusNum);
            myfprintf(verbose, fid, "include: buses %d\n", childBuses);
            myfprintf(verbose, fid, "at indices %d.\n", childBusIndices);
            Aeq(PIdx, indices_P(childBusIndices) ) = -1;   % for P
            Aeq(QIdx, indices_Q(childBusIndices) ) = -1;   % for Q
            for i = 1:length(childBusIndices)
                myfprintf(verbose, fid, "Aeq(%d, P(%d)) = -1\n", PIdx, childBusIndices(i));
                myfprintf(verbose, fid, "Aeq(%d, Q(%d)) = -1\n", QIdx, childBusIndices(i));
            end
        end

        
        % V equations
        % vIdx = parentBusIdx + 2*(N_Area-1);
        vIdx = QIdx + (N_Area-1);
        Aeq( vIdx, indices_v(parentBusIdx) ) = 1;
        myfprintf(verbose, fid, "Aeq(%d, v(%d)) = 1\n", vIdx, parentBusIdx);

        %Return the rows with the list of 'children' buses of 'under' the PARENT of our currentBus:
        %our currentBus will obviously also be included in the list.
        siblingBusesIndices = find(graphDFS_Area_Table.fbus == parentBusNum);
        siblingBuses = graphDFS_Area_Table.tbus(siblingBusesIndices);

        myfprintf(verbose, fid, "The siblings of bus %d\n", currentBusNum);
        myfprintf(verbose, fid, "include these buses: %d\n", siblingBuses)
        myfprintf(verbose, fid, "at indices %d.\n", siblingBusesIndices);
        eldestSiblingIdx = siblingBusesIndices(1);
        eldestSiblingBus = siblingBuses(1);
        myfprintf(verbose, fid,  "which makes bus %d at index %d as the eldest sibling.\n", eldestSiblingBus, eldestSiblingIdx);
        Aeq( vIdx, indices_v_Full( eldestSiblingIdx ) ) = -1;
        myfprintf(verbose, fid, "Aeq(%d, v_Full(%d)) = -1\n", vIdx, eldestSiblingIdx);
        Aeq( vIdx, indices_P(parentBusIdx) ) = 2 * R_Area_Matrix( parentBusNum, currentBusNum );
        myfprintf(verbose, fid, "Aeq(%d, P(%d)) = 2*r(%d, %d).\n", vIdx, parentBusIdx, parentBusNum, currentBusNum);
        Aeq( vIdx, indices_Q(parentBusIdx) ) = 2 * X_Area_Matrix( parentBusNum, currentBusNum );
        myfprintf(verbose, fid, "Aeq(%d, Q(%d)) = 2*x(%d, %d).\n", vIdx, parentBusIdx, parentBusNum, currentBusNum);
        Aeq( vIdx, indices_l(parentBusIdx) ) = ...
            -R_Area_Matrix( parentBusNum, currentBusNum )^2 + ...
            -X_Area_Matrix( parentBusNum, currentBusNum )^2 ;
        myfprintf(verbose, fid, "Aeq(%d, l(%d)) = -r(%d, %d)^2 -x(%d, %d)^2.\n", vIdx, parentBusIdx, parentBusNum, currentBusNum, parentBusNum, currentBusNum);
        
        % beq Formulation
    %     beq = zeros(1, 3*N_Area - 2);
        beq( PIdx ) = ...
            ( 1- 0.5 * CVR_P ) * ...
            ( P_L_Area( currentBusNum ) - P_der_Area( currentBusNum ) );
        myfprintf(verbose, fid, "beq(%d) = (1 - 0.5*CVR_P)*(P_L(%d) - P_der(%d))\n", PIdx, currentBusNum, currentBusNum);
    
        beq( QIdx ) =  ...
            ( 1- 0.5*CVR_Q ) * ...
            ( Q_L_Area( currentBusNum ) - Q_C_Area( currentBusNum ) );
        myfprintf(verbose, fid, "beq(%d) = (1 - 0.5*CVR_Q)*(Q_L(%d) - Q_C(%d))\n", QIdx, currentBusNum, currentBusNum);

    end
    
    % substation voltage equation
    Aeq( 3*(N_Area-1) + 1, indices_v_Full(1) ) = 1;
    beq( 3*(N_Area-1) + 1 ) = v2_parent_Area;
    
    % DER equation addition
    Table_DER = zeros(numBusesWithDERs_Area, 5);
    
    for k22 = 1 : size(busesWithDERs_Area,  1)
        Aeq( indices_Q( Table_Area_Table.tbus == busesWithDERs_Area(k22) ), end+1 ) = 1;
        
        %setting other parameters for DGs:
        Table_DER(k22, 2) = size(Aeq, 2);
        
        % slope kq definiton:
        Table_DER(k22, 3) = 2*ub_Q_onlyDERbuses_Area(k22)/(V_max-V_min); % Qmax at Vmin, and vice versa     
        % Q_ref, V_ref definition:
        Table_DER(k22, 4) = 0.00;  %Qref
        Table_DER(k22, 5) = 1.00;  %Vref
    end
    
    % Table_DER_Table = array2table(Table_DER, 'VariableNames', {'Idx', 'DG_parameter', 'Slope_kq', 'Q_ref', 'V_ref'});
    % mydisplay(verbose, Table_DER_Table)
    
    if fileOpenedFlag
        fclose(fid);
    end
    
    % calling linear solution for intial point
    [V_linear_Area, x_linear_Area, Table_linear_Area, Volttable_linear_Area] = ...
        singlephaselin(v2_parent_Area, S_connection_Area, isLeaf_Area, ...
        Area, N_Area, graphDFS_Area, graphDFS_Area_Table, R_Area_Matrix, X_Area_Matrix, ...
        P_L_Area, Q_L_Area, Q_C_Area, P_der_Area, ...
        busesWithDERs_Area, lb_Q_onlyDERbuses_Area, ub_Q_onlyDERbuses_Area, itr);
    
    Pflow0_Area = x_linear_Area( 1 : N_Area-1 );
    Qflow0_Area = x_linear_Area( N_Area : 2*(N_Area-1) );
    V0_Area =  V_linear_Area;
    Qder0_Area = x_linear_Area( ( end- numBusesWithDERs_Area + 1) : end);
    
    for currentBusNum = 2 : N_Area
        parentBusIdx = find(graphDFS_Area_Table.tbus == currentBusNum);
        siblingBusesIndices = find(parentBusNum == graphDFS_Area_Table.fbus);
        Iflow0( Table_linear_Area(parentBusIdx, 3) )= x_linear_Area(Table_linear_Area(parentBusIdx,3))^2+x_linear_Area(Table_linear_Area(parentBusIdx,4))^2 / x_linear_Area( Volttable_linear_Area( siblingBusesIndices(1) ) );
    end
    
    x0_Area = [Pflow0_Area; Qflow0_Area; Iflow0'; V0_Area; Qder0_Area];
    
    % Definig Limits
    
    lb_Area(1,1) = -1500;                                        % this is to limit the power flow going reverse at the substation
    lb_Area(2:2*(N_Area-1),1)= (-1500*ones(2*(N_Area-1)-1,1));       % P Q limit
    lb_Area(2*(N_Area-1)+1:3*(N_Area-1),1)= zeros((N_Area-1),1);         % I limit
    lb_Area(3*(N_Area-1)+1:4*(N_Area-1)+1,1)= ((V_min^2)*ones(N_Area,1)); % V limit
    
    
    ub_Area(1:2*(N_Area-1),1)= (1500*ones(2*(N_Area-1),1));          % P Q limit
    ub_Area(2*(N_Area-1)+1:3*(N_Area-1),1)= (1500*ones((N_Area-1),1));   % I limit
    ub_Area(3*(N_Area-1)+1:4*(N_Area-1)+1,1)= ((V_max^2)*ones(N_Area,1));      % V limit
    
    lb_Area = [lb_Area; lb_Q_onlyDERbuses_Area];
    ub_Area = [ub_Area; ub_Q_onlyDERbuses_Area];
    
    %  Optimization - 
    options = optimoptions('fmincon', 'Display', 'off', 'MaxFunctionEvaluations', 100000000, 'Algorithm', 'sqp');
    
    % For Loss minimization:
    % mydisplay(verbose, Aeq);
    % mydisplay(verbose, beq);
    
    startSolvingForOptimization = tic;
    
    if itr == 0
        display(graphDFS_Area_Table)
    end

    [x, ~, ~, ~] = fmincon( @(x)objfunTables(x, N_Area, graphDFS_Area_Table.fbus, graphDFS_Area_Table.tbus, indices_l, R_Area_Matrix), ...
                              x0_Area, [], [], Aeq, beq, lb_Area, ub_Area, ...
                              @(x)eqcons(x, Area, N_Area, ...
                              graphDFS_Area_Table.fbus, graphDFS_Area_Table.tbus, indices_P, indices_Q, indices_l, indices_v_Full, ...
                              itr, systemName, numAreas, "verbose", false, "saveToFile", false),...
                              options);
    
    optimizationSolutionTime = toc(startSolvingForOptimization);
    
    time_dist(itr+1, Area) = optimizationSolutionTime;
    
    % Result
    Pall = x(indices_P); %40x1
    Qall = x(indices_Q); %40x1
    Sall = complex(Pall, Qall); %40x1
    P1 = Pall(1); %1x1
    Q1 = Qall(1); %1x1
    S_parent_Area = complex(P1, Q1);  %1x1  % In Pu

    
    Dec_Var_Q = x(end - numBusesWithDERs_Area + 1 : end);

    decisionVars_Area = zeros(N_Area, 1);
    
    for currentBusWithDERIdx = 1 : numBusesWithDERs_Area
        decisionVars_Area( busesWithDERs_Area(currentBusWithDERIdx) ) = Dec_Var_Q(currentBusWithDERIdx);
    end
    
   
    v2_Area = zeros(N_Area, 1);
    S_child_Area = zeros(N_Area, 1);

    for j = 1:size(Table_Area,1)
        v2_Area( Table_Area_Table.tbus(j) ) = x( end - N_Area + 1 - numBusesWithDERs_Area + j);
        S_child_Area( Table_Area_Table.tbus(j) - 1 ) = Sall(j);
    end
    
    v2_Area(1) = v2_parent_Area;
    
    microIterationLosses(itr + 1, Area) = P1 + sum(P_der_Area) - sum(P_L_Area);

end