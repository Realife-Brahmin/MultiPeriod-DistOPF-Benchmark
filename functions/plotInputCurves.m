function plotInputCurves(sysInfo, simInfo, varargin)
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

    f = figure('visible', showPlots);
    yyaxis left;
    plot(t, loadShapePV, '-s', 'Color', [0.85, 0.325, 0.098], 'MarkerFaceColor', [0.85, 0.325, 0.098], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Solar Irradiance');
    hold on;
    plot(t, loadShape, '-^', 'Color', [0.7, 0.7, 0], 'MarkerFaceColor', [0.7, 0.7, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Loading Factor');
    ylabel('Loading Factor [dimensionless]', 'Color', 'k');
    ax = gca;
    ax.YColor = 'k';  % Black for left y-axis

    yyaxis right;
    plot(t, costArray, '-d', 'Color', [0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Cost');
    ylabel('Cost [\$/kWh]', 'Color', [0, 0.5, 0]);
    ax = gca;
    ax.YColor = [0, 0.5, 0];  % Green for right y-axis

    xlabel('Time Period t');
    legend('show', 'Location', 'northwest');
    grid on;
    ax.GridAlpha = 0.5;  % Set grid transparency
    grid minor;
    xticks(1:T);  % Set x-ticks to integers only

    if savePlots
        folderPath = fullfile(wdSim, 'processedData', sysInfo.systemName, 'numAreas', num2str(sysInfo.numAreas));
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        ext = ".png";
        filename = fullfile(folderPath, strcat('InputCurves_Horizon_', num2str(T), ext));
        saveas(f, filename, 'png');
    end

    if showPlots
        set(f, 'Visible', 'on');
    end

    hold off;
end
