function plotInputCurves(sysInfo, simInfo, varargin)
    % plotInputCurves Plots cost, PV load shape, and load shape on a single figure
    % with two different y-axes. Optionally saves and/or shows the plot.

    p = inputParser;
    addParameter(p, 'savePlots', true, @islogical);
    addParameter(p, 'showPlots', false, @islogical);
    parse(p, varargin{:});
    
    costArray = simInfo.costArray;
    loadShapePV = sysInfo.loadShapePV;
    loadShape = sysInfo.loadShape;

    T = simInfo.T;
    savePlots = p.Results.savePlots;
    showPlots = p.Results.showPlots;
    wdSim = simInfo.wdSim;
    t = 1:length(costArray);

    % Adjust visibility based on the 'showPlots' flag
    visibleSetting = 'off';
    if showPlots
        visibleSetting = 'on';
    end

    % Create figure with controlled visibility
    f = figure('Visible', visibleSetting);

    % First y-axis (left side)
    yyaxis left;
    plot(t, loadShapePV, '-s', 'Color', [0.85, 0.325, 0.098], 'MarkerFaceColor', [0.85, 0.325, 0.098], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Solar Irradiance');
    hold on;
    plot(t, loadShape, '-^', 'Color', [0.7, 0.7, 0], 'MarkerFaceColor', [0.7, 0.7, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Loading Factor');
    ylabel('Loading Factor [dimensionless]', 'Color', 'k');

    % Second y-axis (right side)
    yyaxis right;
    plot(t, costArray, '-d', 'Color', [0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Cost');
    ylabel('Cost [\$/kWh]', 'Color', [0, 0.5, 0]);

    % Set y-axis and x-axis limits
    yyaxis left;
    ylim([0, max(loadShapePV)*1.1]);  % Slightly above max for aesthetic spacing
    yyaxis right;
    ylim([0, max(costArray)*1.1]);

    xlim([1, T]);

    % Labels, title, and legend
    xlabel('Time Period t');
    legend('show', 'Location', 'northwest');
    grid minor;

    % Save the figure if requested
    if savePlots
        folderPath = fullfile(wdSim, 'processedData', sysInfo.systemName, 'numAreas', num2str(sysInfo.numAreas));
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        ext = ".png";
        filename = fullfile(folderPath, strcat('InputCurves_Horizon_', num2str(T)));
        saveas(f, filename, 'png');  % Ensure the format is explicitly specified
    end

    % If showing the plot, update the visibility dynamically
    if showPlots
        set(f, 'Visible', 'on');
    end

    % hold off;
end
