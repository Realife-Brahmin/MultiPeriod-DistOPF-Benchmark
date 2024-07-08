function plotInputCurves(sysInfo, simInfo, varargin)
    % plotInputCurves Plots cost, PV load shape, and load shape on a single figure
    % with two different y-axes. Optionally saves and/or shows the plot.

    p = inputParser;
    addParameter(p, 'savePlots', true, @islogical);
    addParameter(p, 'figName', 'plot.png', @ischar);
    addParameter(p, 'showPlots', true, @islogical);
    addParameter(p, 'wdSim', './processedData', @ischar);
    parse(p, varargin{:});
    
    costArray = simInfo.costArray;
    loadShapePV = sysInfo.loadShapePV;
    loadShape = sysInfo.loadShape;
    battstring = simInfo.battstring;
    T = simInfo.T;
    savePlots = p.Results.savePlots;
    figName = strcat("InputCurves_Horizon_", num2str(T), "_", battstring);
    showPlots = p.Results.showPlots;
    wdSim = p.Results.wdSim;

    % t vector
    t = 1:length(costArray);

    % Create figure
    figure('Visible', 'off');

    % First y-axis (left side)
    yyaxis left;
    plot(t, loadShapePV, '-s', 'Color', [0.85, 0.325, 0.098], 'MarkerFaceColor', [0.85, 0.325, 0.098], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Solar Irradiance');
    hold on;
    plot(t, loadShape, '-^', 'Color', [0.7, 0.7, 0], 'MarkerFaceColor', [0.7, 0.7, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Loading Factor');
    ylabel('Loading Factor [dimensionless]', 'Color', 'k');

    % Second y-axis (right side)
    yyaxis right;
    % display(costArray)
    plot(t, costArray, '-d', 'Color', [0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Cost');
    ylabel('Cost [\$/kWh]', 'Color', 'k');  % Use backslash to escape the dollar sign in ylabel

    % Set y-axis limits
    ylim([min([loadShapePV, loadShape], [], 'all')*0.97, max([loadShapePV, loadShape], [], 'all')*1.03]);
    yyaxis right;
    ylim([min(costArray)*0.97, max(costArray)*1.03]);

    % Set x-axis limits
    xlim([1 max(t)*1.03]);

    % Labels, title, and legend
    xlabel('Time Period t');
    % title('t-Series Data Comparison');
    legend('show', 'Location', 'northwest');
    grid on;
    ax = gca;
    ax.GridColor = [0.8, 0.8, 0.8];
    ax.GridAlpha = 0.3;
    ax.MinorGridAlpha = 0.1;  % Lighter minor grid lines

    % Option to save the figure
    if savePlots
        % folderPath = fullfile(wdSim, 'plots');
             
        folderPath = strcat(wdSim, filesep, "processedData", filesep, sysInfo.systemName, ...
            filesep, "numAreas_", num2str(sysInfo.numAreas));
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        print(fullfile(folderPath, figName), '-dpng');
    end

    % Option to show the figure
    if showPlots
        set(gcf, 'Visible', 'on');
    end

    hold off;
end
