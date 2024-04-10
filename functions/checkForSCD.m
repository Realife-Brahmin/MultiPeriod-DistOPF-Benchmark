function checkForSCD(sysInfo, simInfo, areaInfo, T, x, varargin)
    %checkForSCD Plots Charging, Discharging, and SOC for each battery in 
    % each area to check for simultaneous charging and discharging (SCD).

    % Parse input arguments
    p = inputParser;
    % addParameter(p, 'savePlots', false, @(x) islogical(x) || ismember(lower(x), {'true', 'false'}));
    addParameter(p, 'savePlots', false, @(x) islogical(x) || x == 1 || x == 0 || ismember(lower(x), {'true', 'false'}));

    % addParameter(p, 'showPlots', false, @(x) islogical(x) || ismember(lower(x), {'true', 'false'}));
    addParameter(p, 'showPlots', false, @(x) islogical(x) || x == 1 || x == 0 || ismember(lower(x), {'true', 'false'}));

    kVA_B = sysInfo.kVA_B;
    parse(p, varargin{:});
    savePlots = p.Results.savePlots;
    showPlots = p.Results.showPlots;

    nBatt_Area = areaInfo.nBatt_Area;
    fprintf('Checking for SCD for Area %d\n', areaInfo.Area);

    wineRed = [0.7, 0.0, 0.3];
    darkGreen = [0.0, 0.4, 0.0];
    gptPurple = [0.5, 0.2, 0.7];
    
    if showPlots
        visible = 'on';
    else
        visible = 'off';
    end
    
    alpha = simInfo.alpha;
    gamma = simInfo.gamma;
    threshold = 0.001; % Values below this threshold are treated as zero
    
    for batt_num = 1:nBatt_Area
        f = figure('visible', visible);        
        SOC_Max = areaInfo.E_onlyBattBusesMax_Area(batt_num);
        
        % Charging and Discharging subplot
        subplot(2, 1, 1);
        indices_Pc = getIndicesT(areaInfo.indices_Pcj, batt_num);
        indices_Pd = getIndicesT(areaInfo.indices_Pdj, batt_num); 
        
        Pc_1toT_kW = x(indices_Pc)*kVA_B;
        Pd_1toT_kW = -x(indices_Pd)*kVA_B;

        % Apply threshold
        Pc_1toT_kW(abs(Pc_1toT_kW) < threshold) = 0;
        Pd_1toT_kW(abs(Pd_1toT_kW) < threshold) = 0;

        % bar(1:T, x(indices_Pc)*kVA_B, 'FaceColor', darkGreen); 
        bar(1:T, Pc_1toT_kW, 'FaceColor', darkGreen); 
        hold on;
        % bar(1:T, -x(indices_Pd)*kVA_B, 'FaceColor', wineRed);
        bar(1:T, Pd_1toT_kW, 'FaceColor', wineRed);
        hold off;
        
        % keyboard;
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

        % Save plot if requested
        if savePlots
            % folderName = strcat("processedData", filesep, sysInfo.systemName, filesep, "numAreas_", num2str(sysInfo.numAreas), filesep, "area", num2str(areaInfo.Area), filesep, "BatteryVariables");
            battstring = simInfo.battstring;

            folderName = strcat("processedData", filesep, sysInfo.systemName, ...
                filesep, "numAreas_", num2str(sysInfo.numAreas), filesep, ...
                "area", num2str(areaInfo.Area), filesep, "BatteryVariables", filesep, ...
                "Horizon_", num2str(T), filesep, battstring);

            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end

            filename = strcat(folderName, filesep, ...
                "macroItr_", num2str(simInfo.macroItr+1), "_Battery_", num2str(batt_num), ...
                "_alpha_", num2str(alpha), "_gamma_", num2str(gamma),  ".png");
            saveas(f, filename);
        end
    end
end


