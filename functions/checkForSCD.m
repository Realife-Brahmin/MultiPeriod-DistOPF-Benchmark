function checkForSCD(areaInfo, T, x)
    %checkForSCD Plots Charging, Discharging, and SOC for each battery in 
    % each area to check for simultaneous charging and discharging (SCD).

    nBatt_Area = areaInfo.nBatt_Area;
    fprintf('Checking for SCD for Area %d\n', areaInfo.Area);

    for i = 1:nBatt_Area
        figure; % New figure for each battery
        
        % Charging subplot
        subplot(3, 1, 1);
        indices_Pc = getIndicesT(areaInfo.indices_Pcj, i); 
        plot(1:T, x(indices_Pc)* 1000); % convert to kW
        title(sprintf('Area %d Battery %d Charging', areaInfo.Area, i));
        xlabel('Time Interval Number');
        ylabel('kW');
        
        % Discharging subplot
        subplot(3, 1, 2);
        indices_Pd = getIndicesT(areaInfo.indices_Pdj, i; 
        plot(1:T, x(indices_Pd)* 1000);% convert to kW
        title(sprintf('Area %d Battery %d Discharging', areaInfo.Area, i));
        xlabel('Time Interval Number');
        ylabel('kW');
        
        % SOC subplot
        subplot(3, 1, 3);
        indices_B = getIndicesT(areaInfo.indices_Bj, i); 
        plot(1:T, x(indices_B)* 100); % convert to %
        title(sprintf('Area %d Battery %d SOC', areaInfo.Area, i));
        xlabel('Time Interval Number');
        ylabel('%');
    end
    
end
