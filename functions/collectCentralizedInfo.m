function sysInfo = collectCentralizedInfo(sysInfo, simInfo)
    T = simInfo.T;
    N = sysInfo.N;
    m = sysInfo.m;
    kV_B = sysInfo.kV_B;
    kVA_B = sysInfo.kVA_B;
    numAreas = sysInfo.numAreas;
    
    nDER = sysInfo.nDER;
    nBatt = sysInfo.nBatt;
    isRoot = sysInfo.isRoot;

    [sysInfo.bus1] = deal( zeros(N, 1) );
    [sysInfo.fb, sysInfo.tb] = deal( zeros(m, 1) );
    [sysInfo.P_L_1toT, sysInfo.Q_L_1toT] = deal( zeros(m, T) );
    sysInfo.V_1toT = zeros(N, T);
    sysInfo.Q_C_Full = zeros(N, 1);
    [sysInfo.Sder, sysInfo.Pmpp] = deal( zeros(nDER, 1) );

    [sysInfo.pD_1toT, sysInfo.qD_1toT] = deal( zeros(nDER, T) );

    [sysInfo.B0] = deal( zeros(nBatt, 1) );

    [sysInfo.Sbatt, sysInfo.Pbatt] = deal( zeros(nBatt, 1) );
    [sysInfo.B_1toT, sysInfo.Pd_1toT, sysInfo.Pc_1toT, sysInfo.qD_1toT] = deal( zeros(nBatt, T) );

    for areaNum = 1:numAreas
        
        areaInfo = sysInfo.Area{areaNum};
        N_Area = areaInfo.N_Area;
        m_Area = areaInfo.m_Area;
        
        bus1_expanded = areaInfo.bus_Actual; % How to get this?
        fb_expanded = areaInfo.fb_Actual;
        tb_expanded = areaInfo.tb_Actual;
        
        numChildAreas = sysInfo.numChildAreas(areaNum);

        if ~isRoot(areaNum)
            bus1 = bus1_expanded(3:end); % good for copying voltages.
            bus1_Area = 3:N_Area;
            bus1_ownLoads = bus1_expanded(3:end-numChildAreas); % good for copying loads
            bus1_ownLoads_Area = 3:N_Area-numChildAreas;
            fb = fb_expanded(2:end);
            fb_Area = 2:m_Area;
            tb = tb_expanded(2:end);
            tb_Area = 2:m_Area;
        else
            bus1 = bus1_expanded; % good for copying voltages.
            bus1_Area = 1:N_Area; 
            bus1_ownLoads = bus1_expanded(1:end-numChildAreas); % good for copying loads
            bus1_ownLoads_Area = 1:N_Area-numChildAreas;
            fb = fb_expanded;
            fb_Area = 1:m_Area;
            tb = tb_expanded;
            tb_Area = 1:m_Area;
        end
        
        sysInfo.P_L_1toT(bus1_ownLoads, 1:T) = areaInfo.P_L_Area_1toT(bus1_ownLoads_Area, 1:T);

        sysInfo.Q_L_1toT(bus1_ownLoads, 1:T) = areaInfo.Q_L_Area_1toT(bus1_ownLoads_Area, 1:T);

        sysInfo.V_1toT(bus1, 1:T) = areaInfo.V_Area_1toT(bus1_Area, 1:T);

        sysInfo.Q_C_Full(bus1) = areaInfo.Q_C_Area(bus1_Area);

        numDERBus = areaInfo.DERBusNums_Actual;
        busesWithDERs_Area = areaInfo.busesWithDERs_Area;
        sysInfo.Sder(numDERBus) = areaInfo.S_onlyDERbuses_Area;
        sysInfo.Pmpp(numDERBus) = areaInfo.Pmpp_Area(busesWithDERs_Area); % FULL
        sysInfo.pD_1toT(numDERBus, 1:T) = areaInfo.P_der_Area_1toT(busesWithDERs_Area, 1:T); % FULL
        sysInfo.qD_1toT(numDERBus, 1:T) = areaInfo.qD_Area_1toT;

        Batt_percent = simInfo.Batt_percent;
        if Batt_percent > 0
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
    
    end
    
    sysInfo.P_L_Total_1toT = sum(sysInfo.P_L_1toT);
    sysInfo.P_L_Total_allT = sum(sysInfo.P_L_Total_1toT);

    sysInfo.Q_L_Total_1toT = sum(sysInfo.Q_L_1toT);
    sysInfo.Q_L_Total_allT = sum(sysInfo.Q_L_Total_1toT);

    sysInfo.PSubs_1toT = sysInfo.Area{1}.P_Area_1toT(1, 1:T);
    sysInfo.PSubs_allT = sum(sysInfo.PSubs_1toT);
    sysInfo.QSubs_1toT = sysInfo.Area{1}.Q_Area_1toT(1, 1:T);
    sysInfo.QSubs_allT = sum(sysInfo.QSubs_1toT);
    
    if Batt_percent > 0
        sysInfo.Pd_Total_1toT = sum(sysInfo.Pd_1toT);
        sysInfo.Pd_Total_allT = sum(sysInfo.Pd_Total_1toT);
        sysInfo.Pc_Total_1toT = sum(sysInfo.Pc_1toT);
        sysInfo.Pc_Total_allT = sum(sysInfo.Pc_Total_1toT);
        sysInfo.Pdc_Total_1toT = sum(sysInfo.Pd_1toT - sysInfo.Pc_1toT);
        sysInfo.Pdc_Total_allT = sum(sysInfo.Pdc_Total_1toT);
        sysInfo.qB_Total_1toT = sum(sysInfo.qB_1toT);
        sysInfo.qB_Total_allT = sum(sysInfo.qB_Total_1toT);
    else
        sysInfo.Pd_Total_1toT = zeros(T, 1);
        sysInfo.Pd_Total_allT = 0;
        sysInfo.Pc_Total_1toT = zeros(T, 1);
        sysInfo.Pc_Total_allT = 0;
        sysInfo.Pdc_Total_1toT = zeros(T, 1);
        sysInfo.Pdc_Total_allT = 0;
        sysInfo.qB_Total_1toT = zeros(T, 1);
        sysInfo.qB_Total_allT = 0;
    end

    sysInfo.pD_Total_1toT = sum(sysInfo.pD_1toT);
    sysInfo.pD_Total_allT = sum(sysInfo.pD_Total_1toT);
    sysInfo.qD_Total_1toT = sum(sysInfo.qD_1toT);
    sysInfo.qD_Total_allT = sum(sysInfo.qD_Total_1toT);
    
    sysInfo.QC_Total_1toT = repmat(sum(sysInfo.Q_C_Full), 1, T);
    sysInfo.QC_Total_allT = sum(sysInfo.QC_Total_1toT);
    
    sysInfo.p_Total_1toT = sysInfo.pD_Total_1toT + sysInfo.Pdc_Total_1toT;
    sysInfo.p_Total_allT = sysInfo.pD_Total_allT + sysInfo.Pdc_Total_allT;
    
    sysInfo.qGED_Total_1toT = sysInfo.qD_Total_1toT + sysInfo.qB_Total_1toT;
    sysInfo.qGED_Total_allT = sysInfo.qD_Total_allT + sysInfo.qB_Total_allT;
    sysInfo.q_Total_1toT = sysInfo.qGED_Total_1toT + sysInfo.QC_Total_1toT;
    sysInfo.q_Total_allT = sysInfo.qGED_Total_allT + sysInfo.QC_Total_allT;
    
    if Batt_percent > 0
        sysInfo.B_violation = sysInfo.B0 - sysInfo.B_1toT(:, T);
        sysInfo.B_violation_abs = abs(sysInfo.B_violation);
        sysInfo.B_violation_Total = sum(sysInfo.B_violation);
        sysInfo.B_violation_abs_Total = sum(sysInfo.B_violation_abs);
    else
        sysInfo.B_violation = 0;
        sysInfo.B_violation_abs = 0;
        sysInfo.B_violation_Total = 0;
        sysInfo.B_violation_abs_Total = 0;
    end
    
    sysInfo.P_scd_1toT = min(sysInfo.Pd_1toT, sysInfo.Pc_1toT);
    sysInfo.P_scd_Total_1toT = sum(sysInfo.P_scd_1toT);
    sysInfo.P_scd_Total_allT = sum(sysInfo.P_scd_Total_1toT);
    
    sysInfo.P_batt_abs_1toT = max(sysInfo.Pd_1toT, sysInfo.Pc_1toT);
    sysInfo.P_batt_abs_Total_1toT = sum(sysInfo.P_batt_abs_1toT);
    sysInfo.P_batt_abs_Total_allT = sum(sysInfo.P_batt_abs_Total_1toT);

    sysInfo.qB_abs_1toT = abs(sysInfo.qB_1toT);
    sysInfo.qB_abs_Total_1toT = sum(sysInfo.qB_abs_1toT);
    sysInfo.qB_abs_Total_allT = sum(sysInfo.qB_abs_Total_1toT);

    loadShapePV = simInfo.pvCoeffVals;
    sysInfo.loadShapePV = loadShapePV;
    loadShape = simInfo.lambdaVals;
    sysInfo.loadShape = loadShape;

end