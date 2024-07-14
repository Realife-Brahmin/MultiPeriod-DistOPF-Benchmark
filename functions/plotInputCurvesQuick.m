function plotInputCurvesQuick(T, wdSim)
    % Fetch the necessary data
    % wd = pwd();
    [pvCoeffVals, lambdaVals, ~, ~, costArray] = inputForecastData(fullfile(wdSim, "rawData/"), T, 1.0);
    
    % Time vector
    t = 1:T;
    
    % Create figure
    f = figure('visible', 'off');
    
    % First y-axis (left side)
    yyaxis left;
    plot(t, pvCoeffVals, '-s', 'Color', [0.85, 0.325, 0.098], 'MarkerFaceColor', [0.85, 0.325, 0.098], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Solar Irradiance');
    hold on;
    plot(t, lambdaVals, '-^', 'Color', [0.7, 0.7, 0], 'MarkerFaceColor', [0.7, 0.7, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Loading Factor');
    ylabel('Loading Factor [dimensionless]');
    ax = gca;
    ax.YColor = 'k';  % Black for left y-axis

    % Second y-axis (right side)
    yyaxis right;
    plot(t, costArray, '-d', 'Color', [0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Cost');
    ylabel('Cost [\$/kWh]');
    ax = gca;
    ax.YColor = [0, 0.5, 0];  % Green for right y-axis

    xlabel('Time Period t');
    legend('show', 'Location', 'northwest');
    grid minor;
    xticks(1:T);  % Set x-ticks to integers only
    ax.MinorGridLineStyle = '-';
    ax.MinorGridColor = [0.8 0.8 0.8];
    ax.MinorGridAlpha = 0.5;
    ax.YMinorGrid = 'on';
    ax.YMinorTick = 'on';

    % Save the plot in the current directory
    filename = fullfile(wdSim, strcat('InputCurves_Horizon_', num2str(T), '.png'));
    saveas(f, filename, 'png');
    
    hold off;
end
