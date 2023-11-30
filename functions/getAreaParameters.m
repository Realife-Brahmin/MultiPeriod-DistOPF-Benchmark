function areaInfo = getAreaParameters(simInfo, Area, busDataTable_pu_Area, branchDataTable, R_Area, X_Area, varargin)
    %GETAREAPARAMETERS Extracts area-specific power system parameters.
    % 
    % Synopsis:
    %   areaInfo = getAreaParameters(busDataTable_pu_Area, branchDataTable, varargin)
    %
    % Description:
    %   Given the data tables for buses and branches, this function extracts various power 
    %   system parameters specific to an area.
    %
    % Input:
    %   - busDataTable_pu_Area: Table containing bus-specific data.
    %   - branchDataTable: Table containing branch-specific data.
    %   - varargin: Name-Value pairs for 'chargeToPowerRatio', 'soc_min', and 'soc_max'
    %
    % Output:
    %   - areaInfo: Structure containing various power system parameters of the area.
    %
% <<<<<<< HEAD
% 
% =======
% >>>>>>> main
    chargeToPowerRatio = simInfo.chargeToPowerRatio;
    soc_min = simInfo.soc_min;
    soc_max = simInfo.soc_max;
    p = inputParser;
    % addParameter(p, 'chargeToPowerRatio', 4, @isnumeric);
    % addParameter(p, 'soc_min', 0.30, @isnumeric);
    % addParameter(p, 'soc_max', 0.95, @isnumeric);
    % parse(p, varargin{:});
    
    % chargeToPowerRatio = p.Results.chargeToPowerRatio;
    % soc_min = p.Results.soc_min;
    % soc_max = p.Results.soc_max;
    
    areaInfo.Area = Area;
    areaInfo.N_Area = length(busDataTable_pu_Area.bus);
    areaInfo.m_Area = length(branchDataTable.fb);
    areaInfo.fb_Area = branchDataTable.fb;
    areaInfo.tb_Area = branchDataTable.tb;

    % areaInfo.P_L_Area = busDataTable_pu_Area.P_L;
    % keyboard;
    areaInfo.P_L_Area_1toT = busDataTable_pu_Area.P_L_1toT;

    % areaInfo.Q_L_Area = busDataTable_pu_Area.Q_L;
    areaInfo.Q_L_Area_1toT = busDataTable_pu_Area.Q_L_1toT;

    areaInfo.Q_C_Area = busDataTable_pu_Area.Q_C;
    % areaInfo.P_der_Area = busDataTable_pu_Area.P_der;
    areaInfo.P_der_Area0_1toT = busDataTable_pu_Area.P_der0_1toT;
    areaInfo.P_der_Area_1toT = busDataTable_pu_Area.P_der_1toT;
    
    areaInfo.S_der_Area0 = busDataTable_pu_Area.S_der0;
    areaInfo.S_der_Area = busDataTable_pu_Area.S_der;

    areaInfo.S_battMax_Area0 = busDataTable_pu_Area.S_der0;
    areaInfo.S_battMax_Area = busDataTable_pu_Area.S_der;

    P_der0_1toT = busDataTable_pu_Area.P_der0_1toT;
    % areaInfo.P_der0_1toT = P_der0_1toT;
    P_der_1toT = busDataTable_pu_Area.P_der_1toT;
    % areaInfo.P_der_1toT = P_der_1toT;

    areaInfo.P_battMax_Area0 = P_der0_1toT(:, 1);
    areaInfo.P_battMax_Area = P_der_1toT(:, 1);

    areaInfo.Emax_batt_Area0 = chargeToPowerRatio .* areaInfo.P_battMax_Area0;
    areaInfo.Emax_batt_Area = chargeToPowerRatio .* areaInfo.P_battMax_Area;

    % DER Configuration:
    areaInfo.busesWithDERs_Area0 = find(areaInfo.S_der_Area0); % all non-zero element indices
    areaInfo.busesWithDERs_Area = find(areaInfo.S_der_Area); % all non-zero element indices
    areaInfo.nDER_Area0 = length(areaInfo.busesWithDERs_Area0);
    areaInfo.nDER_Area = length(areaInfo.busesWithDERs_Area);
    areaInfo.busesWithBatts_Area0 = find(areaInfo.S_battMax_Area0);
    areaInfo.busesWithBatts_Area = find(areaInfo.S_battMax_Area);
    areaInfo.nBatt_Area0 = length(areaInfo.busesWithBatts_Area0);
    areaInfo.nBatt_Area = length(areaInfo.busesWithBatts_Area);

    areaInfo.S_onlyDERbuses_Area0 = areaInfo.S_der_Area0(areaInfo.busesWithDERs_Area0);
    areaInfo.S_onlyDERbuses_Area = areaInfo.S_der_Area(areaInfo.busesWithDERs_Area);
    % areaInfo.P_onlyDERbuses_Area = areaInfo.P_der_Area(areaInfo.busesWithDERs_Area);
    P_onlyDERbuses_Area0_1toT = areaInfo.P_der_Area0_1toT(areaInfo.busesWithDERs_Area0);
    P_onlyDERbuses_Area_1toT = areaInfo.P_der_Area_1toT(areaInfo.busesWithDERs_Area);
    areaInfo.P_onlyDERbuses_Area0_1toT = P_onlyDERbuses_Area0_1toT;
    areaInfo.P_onlyDERbuses_Area_1toT = P_onlyDERbuses_Area_1toT;
    P_onlyDERbuses_Area0 = P_onlyDERbuses_Area0_1toT(:, 1);
    P_onlyDERbuses_Area = P_onlyDERbuses_Area_1toT(:, 1);
    

    areaInfo.S_onlyBattBusesMax_Area0 = areaInfo.S_battMax_Area0(areaInfo.busesWithBatts_Area0);
    areaInfo.S_onlyBattBusesMax_Area = areaInfo.S_battMax_Area(areaInfo.busesWithBatts_Area);    
    areaInfo.P_onlyBattBusesMax_Area0 = areaInfo.P_battMax_Area0(areaInfo.busesWithBatts_Area0);
    areaInfo.P_onlyBattBusesMax_Area = areaInfo.P_battMax_Area(areaInfo.busesWithBatts_Area);
    areaInfo.E_onlyBattBusesMax_Area0 = areaInfo.Emax_batt_Area0(areaInfo.busesWithBatts_Area0);
    areaInfo.E_onlyBattBusesMax_Area = areaInfo.Emax_batt_Area(areaInfo.busesWithBatts_Area);

    % Set bounds
    areaInfo.lb_Pc_onlyBattBuses_Area0 = zeros(areaInfo.nBatt_Area0, 1);
    areaInfo.lb_Pc_onlyBattBuses_Area = zeros(areaInfo.nBatt_Area, 1);
    areaInfo.ub_Pc_onlyBattBuses_Area0 = areaInfo.P_onlyBattBusesMax_Area0;
    areaInfo.ub_Pc_onlyBattBuses_Area = areaInfo.P_onlyBattBusesMax_Area;
    areaInfo.lb_Pd_onlyBattBuses_Area0 = zeros(areaInfo.nBatt_Area0, 1);
    areaInfo.lb_Pd_onlyBattBuses_Area = zeros(areaInfo.nBatt_Area, 1);
    areaInfo.ub_Pd_onlyBattBuses_Area0 = areaInfo.P_onlyBattBusesMax_Area0;
    areaInfo.ub_Pd_onlyBattBuses_Area = areaInfo.P_onlyBattBusesMax_Area;
    areaInfo.lb_B_onlyBattBuses_Area0 = soc_min.*areaInfo.E_onlyBattBusesMax_Area0;  % soc_min not defined here
    areaInfo.lb_B_onlyBattBuses_Area = soc_min.*areaInfo.E_onlyBattBusesMax_Area;  % soc_min not defined here
    areaInfo.ub_B_onlyBattBuses_Area0 = soc_max.*areaInfo.E_onlyBattBusesMax_Area0;  % soc_max not defined here
    areaInfo.ub_B_onlyBattBuses_Area = soc_max.*areaInfo.E_onlyBattBusesMax_Area;  % soc_max not defined here

    % areaInfo.lb_qD_onlyDERbuses_Area = -sqrt(areaInfo.S_onlyDERbuses_Area.^2 - areaInfo.P_onlyDERbuses_Area.^2);
    areaInfo.lb_qD_onlyDERbuses_Area0 = -sqrt(areaInfo.S_onlyDERbuses_Area0.^2 - P_onlyDERbuses_Area0.^2);
    areaInfo.lb_qD_onlyDERbuses_Area = -sqrt(areaInfo.S_onlyDERbuses_Area.^2 - P_onlyDERbuses_Area.^2);

    % areaInfo.ub_qD_onlyDERbuses_Area = sqrt(areaInfo.S_onlyDERbuses_Area.^2 - areaInfo.P_onlyDERbuses_Area.^2);
    areaInfo.ub_qD_onlyDERbuses_Area0 = sqrt(areaInfo.S_onlyDERbuses_Area0.^2 - P_onlyDERbuses_Area0.^2);
    areaInfo.ub_qD_onlyDERbuses_Area = sqrt(areaInfo.S_onlyDERbuses_Area.^2 - P_onlyDERbuses_Area.^2);

    areaInfo.lb_qB_onlyBattBuses_Area0 = -sqrt(areaInfo.S_onlyBattBusesMax_Area0.^2 - areaInfo.P_onlyBattBusesMax_Area0.^2);
    areaInfo.lb_qB_onlyBattBuses_Area = -sqrt(areaInfo.S_onlyBattBusesMax_Area.^2 - areaInfo.P_onlyBattBusesMax_Area.^2);
    areaInfo.ub_qB_onlyBattBuses_Area0 = sqrt(areaInfo.S_onlyBattBusesMax_Area0.^2 - areaInfo.P_onlyBattBusesMax_Area0.^2);
    areaInfo.ub_qB_onlyBattBuses_Area = sqrt(areaInfo.S_onlyBattBusesMax_Area.^2 - areaInfo.P_onlyBattBusesMax_Area.^2);

    % Initialize Battery SOCs 
    areaInfo.B0Vals_pu_Area0 = mean([areaInfo.lb_B_onlyBattBuses_Area0, areaInfo.ub_B_onlyBattBuses_Area0], 2);
    areaInfo.B0Vals_pu_Area = mean([areaInfo.lb_B_onlyBattBuses_Area, areaInfo.ub_B_onlyBattBuses_Area], 2);

    % The addition of the R_Area_Matrix and X_Area_Matrix:
    [R_Area_Matrix, X_Area_Matrix] = deal(zeros(areaInfo.N_Area, areaInfo.N_Area));
    
    % Matrix form of R and X in terms of graph
    for ij = 1: areaInfo.m_Area
        R_Area_Matrix(areaInfo.fb_Area(ij), areaInfo.tb_Area(ij)) = R_Area(ij);
        R_Area_Matrix(areaInfo.tb_Area(ij), areaInfo.fb_Area(ij)) = R_Area_Matrix(areaInfo.fb_Area(ij), areaInfo.tb_Area(ij));
        X_Area_Matrix(areaInfo.fb_Area(ij), areaInfo.tb_Area(ij)) = X_Area(ij);
        X_Area_Matrix(areaInfo.tb_Area(ij), areaInfo.fb_Area(ij)) = X_Area_Matrix(areaInfo.fb_Area(ij), areaInfo.tb_Area(ij));
    end

    % Assign these matrices to the areaInfo structure
    areaInfo.R_Area_Matrix = R_Area_Matrix;
    areaInfo.X_Area_Matrix = X_Area_Matrix;

end