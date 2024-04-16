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
    batteryTerminalChargeConstraint = simInfo.batteryTerminalChargeConstraint;
    threshold = 0.001; % Values below this threshold are treated as zero
    
    for batt_num = 1:nBatt_Area
        batt_num_Actual = areaInfo.BattBusNums_Actual(batt_num);
        f = figure('visible', visible);        
        SOC_Max = areaInfo.E_onlyBattBusesMax_Area(batt_num);
        
        % Charging and Discharging subplot
        subplot(2, 1, 1);
        indices_Pc = getIndicesT(areaInfo.indices_Pcj, batt_num);
        indices_Pd = getIndicesT(areaInfo.indices_Pdj, batt_num); 
        
        Pc_1toT_kW = x(indices_Pc)*kVA_B;
        Pd_1toT_kW = x(indices_Pd)*kVA_B;

        % Apply threshold
        Pc_1toT_kW(abs(Pc_1toT_kW) < threshold) = 0;
        Pd_1toT_kW(abs(Pd_1toT_kW) < threshold) = 0;

        bar(1:T, Pc_1toT_kW, 'FaceColor', darkGreen); 
        hold on;
        bar(1:T, -Pd_1toT_kW, 'FaceColor', wineRed);
        hold off;
        
        % Modify y-tick labels for Pd to be positive
        ax = gca;

        % Determine dynamic range for y-axis based on data
        maxCharging = max(Pc_1toT_kW);
        maxDischarging = max(abs(Pd_1toT_kW));
        maxYValue = max(1, max(maxCharging, maxDischarging));
        maxYValue = ceil(maxYValue / 5) * 5; % Round up to the nearest 10 for aesthetics
        % keyboard;
        ylim([-maxYValue maxYValue]);
        yticks(linspace(-maxYValue, maxYValue, 6)); % Set y-ticks at every integer within the range

        % Adjust tick labels and grid
        ax = gca;
        ax.YAxis.Exponent = 0;  % Ensure no scientific notation
        ax.YTickLabel = arrayfun(@(v) num2str(v, '%d'), ax.YTick, 'UniformOutput', false);

        % Turn on minor grid lines
        ax.MinorGridLineStyle = '-';
        ax.MinorGridColor = [0.8 0.8 0.8];
        ax.MinorGridAlpha = 0.5;
        ax.YMinorGrid = 'on';
        ax.YMinorTick = 'on';

        % title(sprintf('Area %d Battery %d Charging and Discharging', areaInfo.Area, batt_num));
        title(sprintf('Area %d Battery %d Charging and Discharging', areaInfo.Area, batt_num_Actual));
        xlabel('Time Interval Number');
        ylabel('[kW]');
        legend('Charging', 'Discharging', 'Location', 'southwest');
        grid on;
        grid minor;
        
        % SOC subplot
        subplot(2, 1, 2);
        indices_B = getIndicesT(areaInfo.indices_Bj, batt_num); 
        B0Val = areaInfo.B0Vals_pu_Area(batt_num);
        b = bar([0, 1:T], [B0Val*100/SOC_Max; x(indices_B)*100/SOC_Max], 'FaceColor', gptPurple, 'BarWidth', 0.8); % adjusted bar width
        
        title(sprintf('Area %d Battery %d SOC', areaInfo.Area, batt_num_Actual));
        xlabel('Time Interval Number');
        ylabel('[\%]');
        legend('Battery State of Charge', 'Location', 'southwest');
        grid on;
        grid minor;

        % Save plot if requested
        if savePlots
            battstring = simInfo.battstring;

            folderName = strcat("processedData", filesep, sysInfo.systemName, ...
                filesep, "numAreas_", num2str(sysInfo.numAreas), filesep, ...
                "BatteryVariables", filesep, ...
                "Horizon_", num2str(T), filesep, battstring);

            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end

            if strcmp(batteryTerminalChargeConstraint, "soft")
                filename = strcat(folderName, filesep, ...
                    "macroItr_", num2str(simInfo.macroItr+1), "_Battery_", num2str(batt_num_Actual), ...
                    "_alpha_", num2str(alpha), "_gamma_", num2str(gamma),  ".png");
            elseif strcmp(batteryTerminalChargeConstraint, "hard")
                filename = strcat(folderName, filesep, ...
                "macroItr_", num2str(simInfo.macroItr+1), "_Battery_", num2str(batt_num_Actual), ...
                "_alpha_", num2str(alpha),  ".png");
            else
                error("floc")
            end


            saveas(f, filename);
        end
    end
end


