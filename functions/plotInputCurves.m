function plotInputCurves(costArray, loadShapePV, loadShape, varargin)
    % plotInputCurves Plots cost, PV load shape, and load shape on a single figure
    % with two different y-axes. Optionally saves and/or shows the plot.
    %
    % Inputs:
    %   costArray - Array of cost values in dollars/kWh
    %   loadShapePV - Array of PV load shapes (dimensionless)
    %   loadShape - Array of load shapes (dimensionless)
    %   Optional:
    %     'saveFig' - Boolean to save the figure (default false)
    %     'figName' - Name of the file to save the figure (default 'figure.png')
    %     'showFig' - Boolean to show the figure (default true)

    p = inputParser;
    addParameter(p, 'saveFig', false, @islogical);
    addParameter(p, 'figName', 'figure.png', @ischar);
    addParameter(p, 'showFig', true, @islogical);
    parse(p, varargin{:});

    saveFig = p.Results.saveFig;
    figName = p.Results.figName;
    showFig = p.Results.showFig;

    % Time vector
    time = 1:length(costArray);

    % Create figure
    figure('Visible', 'off');

    % First y-axis (left side)
    yyaxis left;
    plot(time, loadShapePV, '-s', 'Color', [0.85, 0.325, 0.098], 'MarkerFaceColor', [0.85, 0.325, 0.098], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Solar Irradiance');
    hold on;
    plot(time, loadShape, '-^', 'Color', [0.7, 0.7, 0], 'MarkerFaceColor', [0.7, 0.7, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Loading Factor');
    ylabel('Loading Factor [dimensionless]', 'Color', 'k');

    % Second y-axis (right side)
    yyaxis right;
    plot(time, costArray, '-d', 'Color', [0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Cost');
    ylabel('Cost [\$ /kWh]', 'Color', 'k');  % Properly escape the dollar sign

    % Set y-axis limits
    ylim([min([loadShapePV loadShape], [], 'all')*0.97, max([loadShapePV loadShape], [], 'all')*1.03]);
    yyaxis right;
    ylim([min(costArray)*0.97, max(costArray)*1.03]);

    % Set x-axis limits
    xlim([1 max(time)*1.03]);

    % Labels, title, and legend
    xlabel('Time Period t');
    title('Time-Series Data Comparison');
    legend('show', 'Location', 'northwest');
    grid on;
    ax = gca;
    ax.GridColor = [0.8, 0.8, 0.8];
    ax.GridAlpha = 0.3;
    ax.MinorGridAlpha = 0.1;  % Lighter minor grid lines

    % Option to save the figure
    if saveFig
        print(figName, '-dpng');
    end

    % Option to show the figure
    if showFig
        set(gcf, 'Visible', 'on');
    end

    hold off;
end
