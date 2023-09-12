function [N_Area, m_Area, fb_Area, tb_Area, P_der_Area, P_L_Area, Q_L_Area, Q_C_Area, S_der_Area, S_battMax_Area, P_battMax_Area, Emax_batt_Area, busesWithDERs_Area, nDER_Area, busesWithBatts_Area, nBatt_Area, B0Vals_pu_Area] = getAreaParameters(busDataTable, branchDataTable, lambdaVals, pvCoeffVals, B0Vals_pu_Area, t, varargin)
    %GETAREAPARAMETERS Extract time-variant area-specific power system parameters.
    % 
    % Synopsis:
    %   [N_Area, m_Area, fb_Area, tb_Area, P_der_Area, P_L_Area, Q_L_Area, Q_C_Area, S_der_Area, S_battMax_Area, P_battMax_Area, Emax_batt_Area, busesWithDERs_Area, nDER_Area, busesWithBatts_Area, nBatt_Area, B0Vals_pu_Area] 
    %   = getAreaParameters(busDataTable, branchDataTable, lambdaVals, pvCoeffVals, B0Vals_pu, t, varargin)
    %
    % Description:
    %   Given the data tables for buses and branches, as well as time-variant lambdaVals and 
    %   pvCoeffVals, this function extracts various power system parameters specific to an area 
    %   at the time index 't'.
    %
    % Input:
    %   - busDataTable: Table containing bus-specific data.
    %   - branchDataTable: Table containing branch-specific data.
    %   - lambdaVals: Array of scalar multipliers for loads.
    %   - pvCoeffVals: Array of scalar multipliers for PV power generation.
    %   - B0Vals_pu: Initial battery values (states of charge) as an array.
    %   - t: Time index to extract specific lambda and pvCoeff values.
    %   - varargin: Name-Value pair for 'chargeToPowerRatio'
    %
    % Output:
    %   [see above output descriptions]
    %
    
    p = inputParser;
    addParameter(p, 'chargeToPowerRatio', 4, @isnumeric);
    parse(p, varargin{:});
    
    chargeToPowerRatio = p.Results.chargeToPowerRatio;

    N_Area = length(busDataTable.bus);
    m_Area = length(branchDataTable.fb);
    fb_Area = branchDataTable.fb;
    tb_Area = branchDataTable.tb;

    P_L_Area = lambdaVals(t) * busDataTable.P_L;
    Q_L_Area = lambdaVals(t) * busDataTable.Q_L;
    Q_C_Area = busDataTable.Q_C;
    P_der_Area = pvCoeffVals(t) * busDataTable.P_der;
    S_der_Area = busDataTable.S_der;
    S_battMax_Area = busDataTable.S_der;
    P_battMax_Area = busDataTable.P_der;
    Emax_batt_Area = chargeToPowerRatio .* P_battMax_Area;

    % DER Configuration:
    busesWithDERs_Area = find(S_der_Area); % all non-zero element indices
    nDER_Area = length(busesWithDERs_Area);
    busesWithBatts_Area = find(S_battMax_Area);
    nBatt_Area = length(busesWithBatts_Area);
    B0Vals_pu_Area = reshape(B0Vals_pu_Area(1:nBatt_Area), nBatt_Area, 1);
end
