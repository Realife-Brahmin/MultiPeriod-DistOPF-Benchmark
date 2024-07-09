function plotInputCurves(sysInfo, simInfo, varargin)
    % plotInputCurves Plots cost, PV load shape, and load shape on a single figure
    % with two different y-axes. Optionally saves and/or shows the plot.

    p = inputParser;
    addParameter(p, 'savePlots', true, @islogical);
    addParameter(p, 'showPlots', true, @islogical);
    parse(p, varargin{:});
    
    costArray = simInfo.costArray;
    loadShapePV = sysInfo.loadShapePV;
    loadShape = sysInfo.loadShape;

    T = simInfo.T;
    savePlots = p.Results.savePlots;
    figName = strcat("InputCurves_Horizon_", num2str(T));
    ext = ".png";
    showPlots = p.Results.showPlots;
    wdSim = simInfo.wdSim;
    t = 1:length(costArray);

    % Create figure
    f = figure('visible', showPlots);        

    % First y-axis (left side)
    yyaxis left;
    plot(t, loadShapePV, '-s', 'Color', [0.85, 0.325, 0.098], 'MarkerFaceColor', [0.85, 0.325, 0.098], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Solar Irradiance');
    hold on;
    plot(t, loadShape, '-^', 'Color', [0.7, 0.7, 0], 'MarkerFaceColor', [0.7, 0.7, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Loading Factor');
    ylabel('Loading Factor [dimensionless]', 'Color', 'k');

    % Second y-axis (right side)
    yyaxis right;
    plot(t, costArray, '-d', 'Color', [0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Cost');
    ylabel('Cost [\$/kWh]', 'Color', 'k');  % Use backslash to escape the dollar sign in ylabel

    % Set y-axis limits to be tighter
    yyaxis left;
    ylim([0, 1]);  % Adjust these limits as necessary

    yyaxis right;
    ylim([0, 0.35]);  % Adjust these limits as necessary

    % Set x-axis limits
    xlim([1 max(t)]);

    % Labels, title, and legend
    xlabel('Time Period t');
    legend('show', 'Location', 'northwest');
    grid on;
    ax = gca;
    ax.GridColor = [0.8, 0.8, 0.8];
    ax.GridAlpha = 0.3;
    ax.MinorGridAlpha = 0.1;  % Lighter minor grid lines

    % Option to save the figure
    if savePlots             
        folderPath = strcat(wdSim, filesep, "processedData", filesep, sysInfo.systemName, ...
            filesep, "numAreas_", num2str(sysInfo.numAreas));
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        filename = fullfile(folderPath, strcat(figName, ext));
        saveas(f, filename);
    end

    hold off;
end
