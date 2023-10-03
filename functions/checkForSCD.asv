function checkForSCD(areaInfo, T, x)
    %checkForSCD Plots Charging, Discharging, and SOC for each battery in 
    % each area to check for simultaneous charging and discharging (SCD).

    nBatt_Area = areaInfo.nBatt_Area;
    fprintf('Checking for SCD for Area %d\n', areaInfo.Area);

    wineRed = [0.7, 0.0, 0.3];
    darkGreen = [0.0, 0.4, 0.0];
    gptPurple = [0.5, 0.2, 0.7];

    for batt_num = 1:nBatt_Area
        figure; % New figure for each battery
        
        SOC_Max = areaInfo.E_onlyBattBusesMax_Area(batt_num);
        
        % Charging and Discharging subplot
        subplot(2, 1, 1);
        indices_Pc = getIndicesT(areaInfo.indices_Pcj, batt_num);
        indices_Pd = getIndicesT(areaInfo.indices_Pdj, batt_num); 
        
        bar(1:T, x(indices_Pc)*1000, 'FaceColor', darkGreen); 
        hold on;
        bar(1:T, -x(indices_Pd)*1000, 'FaceColor', wineRed);
        hold off;
        
        % Modify y-tick labels for Pd to be positive
        ax = gca;
        ax.YTickLabel = cellfun(@(v) num2str(abs(str2double(v))), ax.YTickLabel, 'UniformOutput', false);
        
        title(sprintf('Area %d Battery %d Charging and Discharging', areaInfo.Area, batt_num));
        xlabel('Time Interval Number');
        ylabel('[kW]');
        legend('Charging', 'Discharging', 'Location', 'NorthWest');
        grid on;
        grid minor;
        
        % SOC subplot
        subplot(2, 1, 2);
        indices_B = getIndicesT(areaInfo.indices_Bj, batt_num); 
        b = bar(1:T, x(indices_B)*100/SOC_Max, 'FaceColor', gptPurple, 'BarWidth', 0.8); % adjusted bar width
        
        title(sprintf('Area %d Battery %d SOC', areaInfo.Area, batt_num));
        xlabel('Time Interval Number');
        ylabel('[\%]');
        legend('Battery State of Charge', 'Location', 'NorthWest');
        grid on;
        grid minor;
    end
end


