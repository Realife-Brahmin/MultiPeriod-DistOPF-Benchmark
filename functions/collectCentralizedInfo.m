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
        
    end
end