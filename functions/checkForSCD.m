function checkForSCD(areaInfo, T, x)
    %checkForSCD Plots Charging, Discharging, and SOC for each battery in 
    % each area to check for simultaneous charging and discharging (SCD).

    nBatt_Area = areaInfo.nBatt_Area;
    fprintf('Checking for SCD for Area %d\n', areaInfo.Area);

    for batt_num = 1:nBatt_Area
        figure; % New figure for each battery
        
        SOC_Max = areaInfo.E_onlyBattBusesMax_Area(batt_num);
        % Charging subplot
        subplot(3, 1, 1);
        indices_Pc = getIndicesT(areaInfo.indices_Pcj, batt_num); 
        plot(1:T, x(indices_Pc)); % convert to kW
        title(sprintf('Area %d Battery %d Charging', areaInfo.Area, batt_num));
        xlabel('Time Interval Number');
        ylabel('[kW]');
        
        % Discharging subplot
        subplot(3, 1, 2);
        indices_Pd = getIndicesT(areaInfo.indices_Pdj, batt_num); 
        plot(1:T, x(indices_Pd)* 1000);% convert to kW
        title(sprintf('Area %d Battery %d Discharging', areaInfo.Area, batt_num));
        xlabel('Time Interval Number');
        ylabel('[kW]');
        
        % SOC subplot
        subplot(3, 1, 3);
        indices_B = getIndicesT(areaInfo.indices_Bj, batt_num); 
        plot(1:T, x(indices_B)*100/SOC_Max); % convert to %
        title(sprintf('Area %d Battery %d SOC', areaInfo.Area, batt_num));
        xlabel('Time Interval Number');
        ylabel('[\%]');
    end
    
end
