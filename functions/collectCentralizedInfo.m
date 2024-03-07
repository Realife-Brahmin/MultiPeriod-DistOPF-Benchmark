function sysInfo = collectCentralizedInfo(sysInfo, simInfo)
    T = simInfo.T;
    N = sysInfo.N;
    m = sysInfo.m;
    kV_B = sysInfo.kV_B;
    kVA_B = sysInfo.kVA_B;
    numAreas = sysInfo.numAreas;
    
    % sysInfo already has nDER and nBatt
    nDER = sysInfo.nDER;
    nBatt = sysInfo.nBatt;

    [sysInfo.P_L_1toT, sysInfo.Q_L_1toT] = deal( zeros(m, T) );
    sysInfo.V_1toT = zeros(N, T);
    sysInfo.Q_C_Full = zeros(N, 1);
    [sysInfo.busesWithDERs, sysInfo.Sder, sysInfo.Pmpp] = deal( zeros(nDER, 1) );
    [sysInfo.pD_1toT, sysInfo.qD_1toT] = deal( zeros(nDER, T) );

    [sysInfo.busesWithBatts, sysInfo.B0] = deal( zeros(nBatt, 1) );
    [sysInfo.Sbatt, sysInfo.Pbatt] = deal( zeros(nBatt, 1) );
    [sysInfo.B_1toT, sysInfo.Pd_1toT, sysInfo.Pc_1toT, sysInfo.qD_1toT] = deal( zeros(nBatt, T) );

    for areaNum = 1:numAreas
    % for areaNum = [2 3 1 4]

        areaInfo = sysInfo.Area{areaNum};
        
        % error("Have you inserted actual bus1, actual fb and actual tb values for the area?")
        
        % busData, branchData
        bus1 = areaInfo.bus_Actual; % How to get this?
        fb = areaInfo.fb_Actual;
        tb = areaInfo.tb_Actual;
        % sysInfo.PL0(bus1) = areaInfo.
        sysInfo.P_L_1toT(bus1, 1:T) = areaInfo.P_L_Area_1toT;
        sysInfo.Q_L_1toT(bus1, 1:T) = areaInfo.Q_L_Area_1toT;
        % keyboard;
        sysInfo.V_1toT(bus1, 1:T) = areaInfo.V_Area_1toT;
        sysInfo.Q_C_Full(bus1) = areaInfo.Q_C_Area;

        % figure out how to create sysInfo.pD_1toT with only nDER buses

        % figure out how to create busesWithDERSonly for sysInfo
        
        % DER Parameters
        % numDERBus will be like [1 7 22 85] which will only contain values
        % between 1:n_DER_System, so no actual bus numbers
        
        % sysInfo.busesWithDERs = 
        % sysInfo.

        % % keyboard;
        numDERBus = areaInfo.DERBusNums_Actual;
        % areaInfo.S_der_Area % FULL
        busesWithDERs_Area = areaInfo.busesWithDERs_Area;
        sysInfo.Sder(numDERBus) = areaInfo.S_onlyDERbuses_Area;
        sysInfo.Pmpp(numDERBus) = areaInfo.Pmpp_Area(busesWithDERs_Area); % FULL
        sysInfo.pD_1toT(numDERBus, 1:T) = areaInfo.P_der_Area_1toT(busesWithDERs_Area); % FULL
        sysInfo.qD_1toT(numDERBus, 1:T) = areaInfo.qD_Area_1toT;
        
        % keyboard;

        % Battery Parameters
        % numBattBus will be like [1 7 22 85] which will only contain values
        % between 1:n_DER_System, so no actual bus numbers
        % sysInfo.busesWithBatts(numBattBus) = areaInfo.busesWithBatts_Area;
        numBattBus = areaInfo.BattBusNums_Actual;
        busesWithBatts_Area = areaInfo.busesWithBatts_Area;
        sysInfo.Sbatt(numBattBus) = areaInfo.S_battMax_Area(busesWithBatts_Area); % currently named as S_battRated in vald
        sysInfo.Pbatt(numBattBus) = areaInfo.P_battMax_Area(busesWithBatts_Area); % currently named as P_battRated in vald
        sysInfo.B0(numBattBus) = areaInfo.B0Vals_pu_Area; 
        sysInfo.B_1toT(numBattBus, 1:T) = areaInfo.B_Area_1toT;
        sysInfo.Pd_1toT(numBattBus, 1:T) = areaInfo.Pd_Area_1toT;
        sysInfo.Pc_1toT(numBattBus, 1:T) = areaInfo.Pc_Area_1toT;
        sysInfo.qB_1toT(numBattBus, 1:T) = areaInfo.qB_Area_1toT;
    
    end
    
    loadShapePV = simInfo.pvCoeffVals;
    sysInfo.loadShapePV = loadShapePV;
    loadShape = simInfo.lambdaVals;
    sysInfo.loadShape = loadShape;

end