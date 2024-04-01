function sysInfo = getSysParameters(sysInfo, simInfo, busDataTable_pu_Area, branchDataTable, R_Area, X_Area, varargin)
    %GETAREAPARAMETERS Extracts area-specific power system parameters.
    % 
    % Synopsis:
    %   sysInfo = getAreaParameters(busDataTable_pu_Area, branchDataTable, varargin)
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
    %   - sysInfo: Structure containing various power system parameters of the area.
    %

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
    
    % sysInfo.Area = 1;
    sysInfo.N_Area = length(busDataTable_pu_Area.bus);
    sysInfo.m_Area = length(branchDataTable.fb);
    sysInfo.bus_Area = busDataTable_pu_Area.bus;
    sysInfo.bus_Actual = busDataTable_pu_Area.busActual; 
    sysInfo.fb_Area = branchDataTable.fb;
    sysInfo.tb_Area = branchDataTable.tb;
    sysInfo.fb_Actual = branchDataTable.fbActual;
    sysInfo.tb_Actual = branchDataTable.tbActual;
    % sysInfo.P_L_Area = busDataTable_pu_Area.P_L;
    % keyboard;
    sysInfo.P_L_Area_1toT = busDataTable_pu_Area.P_L_1toT;

    % sysInfo.Q_L_Area = busDataTable_pu_Area.Q_L;
    sysInfo.Q_L_Area_1toT = busDataTable_pu_Area.Q_L_1toT;

    sysInfo.Q_C_Area = busDataTable_pu_Area.Q_C;
    % sysInfo.P_der_Area = busDataTable_pu_Area.P_der;
    
    sysInfo.Pmpp_Area0 = busDataTable_pu_Area.Pmpp_Area0;
    sysInfo.Pmpp_Area = busDataTable_pu_Area.Pmpp_Area;
    
    sysInfo.P_der_Area0_1toT = busDataTable_pu_Area.P_der0_1toT;
    sysInfo.P_der_Area_1toT = busDataTable_pu_Area.P_der_1toT;
    
    sysInfo.S_der_Area0 = busDataTable_pu_Area.S_der0;
    sysInfo.S_der_Area = busDataTable_pu_Area.S_der;

    % sysInfo.S_battMax_Area0 = busDataTable_pu_Area.S_der0;
    sysInfo.S_battMax_Area0 = busDataTable_pu_Area.S_batt0;

    % sysInfo.S_battMax_Area = busDataTable_pu_Area.S_der;
    sysInfo.S_battMax_Area = busDataTable_pu_Area.S_batt;


    % P_der0_1toT = busDataTable_pu_Area.P_der0_1toT;
    P_batt0 = busDataTable_pu_Area.P_batt0;
    % sysInfo.P_der0_1toT = P_der0_1toT;
    % P_der_1toT = busDataTable_pu_Area.P_der_1toT;
    P_batt = busDataTable_pu_Area.P_batt;
    % sysInfo.P_der_1toT = P_der_1toT;

    % sysInfo.P_battMax_Area0 = P_der0_1toT(:, 1);
    sysInfo.P_battMax_Area0 = P_batt0;

    % sysInfo.P_battMax_Area = P_der_1toT(:, 1);
    sysInfo.P_battMax_Area = P_batt;


    sysInfo.Emax_batt_Area0 = chargeToPowerRatio .* sysInfo.P_battMax_Area0;
    sysInfo.Emax_batt_Area = chargeToPowerRatio .* sysInfo.P_battMax_Area;

    % DER Configuration:
    sysInfo.busesWithDERs_Area0 = find(sysInfo.S_der_Area0); % all non-zero element indices
    sysInfo.busesWithDERs_Actual0 = sysInfo.bus_Actual(sysInfo.busesWithDERs_Area0);
    sysInfo.busesWithDERs_Area = find(sysInfo.S_der_Area); % all non-zero element indices
    sysInfo.busesWithDERs_Actual = sysInfo.bus_Actual(sysInfo.busesWithDERs_Area);

    sysInfo.nDER_Area0 = length(sysInfo.busesWithDERs_Area0);

    sysInfo.nDER_Area = length(sysInfo.busesWithDERs_Area);
    sysInfo.busesWithBatts_Area0 = find(sysInfo.S_battMax_Area0);
    sysInfo.busesWithBatts_Actual0 = sysInfo.bus_Actual(sysInfo.busesWithBatts_Area0);

    sysInfo.busesWithBatts_Area = find(sysInfo.S_battMax_Area);
    sysInfo.busesWithBatts_Actual = sysInfo.bus_Actual(sysInfo.busesWithBatts_Area);
    % keyboard;

    sysInfo.nBatt_Area0 = length(sysInfo.busesWithBatts_Area0);
    sysInfo.nBatt_Area = length(sysInfo.busesWithBatts_Area);

    sysInfo.S_onlyDERbuses_Area0 = sysInfo.S_der_Area0(sysInfo.busesWithDERs_Area0);
    sysInfo.S_onlyDERbuses_Area = sysInfo.S_der_Area(sysInfo.busesWithDERs_Area);
    % sysInfo.P_onlyDERbuses_Area = sysInfo.P_der_Area(sysInfo.busesWithDERs_Area);
    P_onlyDERbuses_Area0_1toT = sysInfo.P_der_Area0_1toT(sysInfo.busesWithDERs_Area0);
    P_onlyDERbuses_Area_1toT = sysInfo.P_der_Area_1toT(sysInfo.busesWithDERs_Area);
    sysInfo.P_onlyDERbuses_Area0_1toT = P_onlyDERbuses_Area0_1toT;
    sysInfo.P_onlyDERbuses_Area_1toT = P_onlyDERbuses_Area_1toT;
    % P_onlyDERbuses_Area0 = P_onlyDERbuses_Area0_1toT(:, 1);
    busesWithDERs_Area0 = sysInfo.busesWithDERs_Area0;
    busesWithDERs_Area = sysInfo.busesWithDERs_Area;
    P_onlyDERbuses_Area0 = sysInfo.Pmpp_Area0(busesWithDERs_Area0);
    % P_onlyDERbuses_Area = P_onlyDERbuses_Area_1toT(:, 1);
    P_onlyDERbuses_Area = sysInfo.Pmpp_Area(busesWithDERs_Area);


    sysInfo.S_onlyBattBusesMax_Area0 = sysInfo.S_battMax_Area0(sysInfo.busesWithBatts_Area0);
    sysInfo.S_onlyBattBusesMax_Area = sysInfo.S_battMax_Area(sysInfo.busesWithBatts_Area);    
    sysInfo.P_onlyBattBusesMax_Area0 = sysInfo.P_battMax_Area0(sysInfo.busesWithBatts_Area0);
    sysInfo.P_onlyBattBusesMax_Area = sysInfo.P_battMax_Area(sysInfo.busesWithBatts_Area);
    sysInfo.E_onlyBattBusesMax_Area0 = sysInfo.Emax_batt_Area0(sysInfo.busesWithBatts_Area0);
    sysInfo.E_onlyBattBusesMax_Area = sysInfo.Emax_batt_Area(sysInfo.busesWithBatts_Area);

    % Set bounds
    sysInfo.lb_Pc_onlyBattBuses_Area0 = zeros(sysInfo.nBatt_Area0, 1);
    sysInfo.lb_Pc_onlyBattBuses_Area = zeros(sysInfo.nBatt_Area, 1);
    sysInfo.ub_Pc_onlyBattBuses_Area0 = sysInfo.P_onlyBattBusesMax_Area0;
    sysInfo.ub_Pc_onlyBattBuses_Area = sysInfo.P_onlyBattBusesMax_Area;
    sysInfo.lb_Pd_onlyBattBuses_Area0 = zeros(sysInfo.nBatt_Area0, 1);
    sysInfo.lb_Pd_onlyBattBuses_Area = zeros(sysInfo.nBatt_Area, 1);
    sysInfo.ub_Pd_onlyBattBuses_Area0 = sysInfo.P_onlyBattBusesMax_Area0;
    sysInfo.ub_Pd_onlyBattBuses_Area = sysInfo.P_onlyBattBusesMax_Area;
    sysInfo.lb_B_onlyBattBuses_Area0 = soc_min.*sysInfo.E_onlyBattBusesMax_Area0;  % soc_min not defined here
    sysInfo.lb_B_onlyBattBuses_Area = soc_min.*sysInfo.E_onlyBattBusesMax_Area;  % soc_min not defined here
    sysInfo.ub_B_onlyBattBuses_Area0 = soc_max.*sysInfo.E_onlyBattBusesMax_Area0;  % soc_max not defined here
    sysInfo.ub_B_onlyBattBuses_Area = soc_max.*sysInfo.E_onlyBattBusesMax_Area;  % soc_max not defined here

    % sysInfo.lb_qD_onlyDERbuses_Area = -sqrt(sysInfo.S_onlyDERbuses_Area.^2 - sysInfo.P_onlyDERbuses_Area.^2);
    % keyboard;
    sysInfo.lb_qD_onlyDERbuses_Area0 = -sqrt(sysInfo.S_onlyDERbuses_Area0.^2 - P_onlyDERbuses_Area0.^2);
    % sysInfo.lb_qD_onlyDERbuses_Area0 = -sqrt(sysInfo.S_onlyDERbuses_Area0.^2 - sysInfo.Pmpp_Area0.^2);
    sysInfo.lb_qD_onlyDERbuses_Area = -sqrt(sysInfo.S_onlyDERbuses_Area.^2 - P_onlyDERbuses_Area.^2);
    % sysInfo.lb_qD_onlyDERbuses_Area = -sqrt(sysInfo.S_onlyDERbuses_Area.^2 - sysInfo.Pmpp_Area.^2);


    % sysInfo.ub_qD_onlyDERbuses_Area = sqrt(sysInfo.S_onlyDERbuses_Area.^2 - sysInfo.P_onlyDERbuses_Area.^2);
    sysInfo.ub_qD_onlyDERbuses_Area0 = sqrt(sysInfo.S_onlyDERbuses_Area0.^2 - P_onlyDERbuses_Area0.^2);
    % sysInfo.ub_qD_onlyDERbuses_Area0 = sqrt(sysInfo.S_onlyDERbuses_Area0.^2 - sysInfo.Pmpp_Area0.^2);

    sysInfo.ub_qD_onlyDERbuses_Area = sqrt(sysInfo.S_onlyDERbuses_Area.^2 - P_onlyDERbuses_Area.^2);
    % sysInfo.ub_qD_onlyDERbuses_Area = sqrt(sysInfo.S_onlyDERbuses_Area.^2 - sysInfo.Pmpp_Area.^2);


    sysInfo.lb_qB_onlyBattBuses_Area0 = -sqrt(sysInfo.S_onlyBattBusesMax_Area0.^2 - sysInfo.P_onlyBattBusesMax_Area0.^2);
    sysInfo.lb_qB_onlyBattBuses_Area = -sqrt(sysInfo.S_onlyBattBusesMax_Area.^2 - sysInfo.P_onlyBattBusesMax_Area.^2);
    sysInfo.ub_qB_onlyBattBuses_Area0 = sqrt(sysInfo.S_onlyBattBusesMax_Area0.^2 - sysInfo.P_onlyBattBusesMax_Area0.^2);
    sysInfo.ub_qB_onlyBattBuses_Area = sqrt(sysInfo.S_onlyBattBusesMax_Area.^2 - sysInfo.P_onlyBattBusesMax_Area.^2);

    % Initialize Battery SOCs 
    sysInfo.B0Vals_pu_Area0 = mean([sysInfo.lb_B_onlyBattBuses_Area0, sysInfo.ub_B_onlyBattBuses_Area0], 2);
    sysInfo.B0Vals_pu_Area = mean([sysInfo.lb_B_onlyBattBuses_Area, sysInfo.ub_B_onlyBattBuses_Area], 2);

    % The addition of the R_Area_Matrix and X_Area_Matrix:
    [R_Area_Matrix, X_Area_Matrix] = deal(zeros(sysInfo.N_Area, sysInfo.N_Area));
    
    % Matrix form of R and X in terms of graph
    for ij = 1: sysInfo.m_Area
        R_Area_Matrix(sysInfo.fb_Area(ij), sysInfo.tb_Area(ij)) = R_Area(ij);
        R_Area_Matrix(sysInfo.tb_Area(ij), sysInfo.fb_Area(ij)) = R_Area_Matrix(sysInfo.fb_Area(ij), sysInfo.tb_Area(ij));
        X_Area_Matrix(sysInfo.fb_Area(ij), sysInfo.tb_Area(ij)) = X_Area(ij);
        X_Area_Matrix(sysInfo.tb_Area(ij), sysInfo.fb_Area(ij)) = X_Area_Matrix(sysInfo.fb_Area(ij), sysInfo.tb_Area(ij));
    end

    % Assign these matrices to the sysInfo structure
    sysInfo.R_Area_Matrix = R_Area_Matrix;
    sysInfo.X_Area_Matrix = X_Area_Matrix;

end