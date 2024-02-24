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
        sysInfo.Pmpp_Full(bus1) = areaInfo.Pmpp_Area

        sysInfo.pD_Full_1toT(bus1, 1:T) = areaInfo.P_der_Area_1toT;
        % figure out how to create sysInfo.pD_1toT with only nDER buses

        % figure out how to create busesWithDERSonly for sysInfo

    end
end