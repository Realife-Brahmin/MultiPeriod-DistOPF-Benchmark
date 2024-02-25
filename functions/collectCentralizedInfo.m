function sysInfo = collectCentralizedInfo(sysInfo, simInfo)
    T = simInfo.T;
    N = sysInfo.N;
    m = sysInfo.m;
    numAreas = sysInfo.numAreas;
    for areaNum = 1:numAreas
        areaInfo = sysInfo.Area{areaNum};

        error("Have you inserted actual bus1, actual fb and actual tb values for the area?")
        
        bus1 = areaInfo.bus1_Actual;
        % fb_Area = areaInfo.fb_Area;
        fb = areaInfo.fb_Actual;
        % tb_Area = areaInfo.tb_Area;
        tb = areaInfo.tb_Actual;
        % P_L_1toT = sysInfo.P_L_1toT;
        sysInfo.P_L_1toT(bus1, 1:T) = areaInfo.P_L_Area_1toT;
        sysInfo.Q_L_1toT(bus1, 1:T) = areaInfo.Q_L_Area_1toT;
        
        sysInfo.Q_C_Full(bus1) = areaInfo.Q_C_Area
        % sysInfo.Pmpp_Full(bus1) = areaInfo.Pmpp_Area

        % sysInfo.pD_Full_1toT(bus1, 1:T) = areaInfo.P_der_Area_1toT;
        % figure out how to create sysInfo.pD_1toT with only nDER buses

        % figure out how to create busesWithDERSonly for sysInfo
        
        % bus1_DER = areaInfo.bus1_Actual_DER; % this gives me the actual bus number where the DER is placed.
        % This is not useful for non-Full variables for sysInfo for me.
        % I somehow need to get the batt_num_sys values
        
        % numDERBus will be like [1 7 22 85] which will only contain values
        % between 1:n_DER_System, so no actual bus numbers
        sysInfo.Sder(numDERBus) = areaInfo.S_der_Area;
        sysInfo.Pmpp(numDERBus) = areaInfo.Pmpp_Area;
        sysInfo.pD_1toT(numDERBus, 1:T) = areaInfo.P_der_Area_1toT;
        sysInfo.qD_1toT(numDERBus, 1:T) = areaInfo.qD_Area_1toT;
        
        % numBattBus will be like [1 7 22 85] which will only contain values
        % between 1:n_DER_System, so no actual bus numbers
        sysInfo.Sbatt(numBattBus) = areaInfo.S_battMax_Area;
        sysInfo.B0(numBattBus) = areaInfo.B0Vals_pu_Area; 
        sysInfo.B_1toT(numBattBus, 1:T) = areaInfo.B_Area_1toT;
        sysInfo.Pd_1toT(numBattBus, 1:T) = areaInfo.Pd_Area_1toT;
        sysInfo.Pc_1toT(numBattBus, 1:T) = areaInfo.Pc_Area_1toT;
        sysInfo.qB_1toT(numBattBus, 1:T) = areaInfo.qB_Area_1toT;
    
        

    end
    
    % loadShapePV = simInfo.pvCoeffVals;

end