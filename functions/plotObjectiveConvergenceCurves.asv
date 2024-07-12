function plotObjectiveConvergenceCurves(sysInfo, simInfo, varargin)
    % plotObjectiveConvergenceCurves Plots convergence curves for Substation power and cost on a single plot with dual y-axes.
    % Inputs:
    %   sysInfo - System information structure with fields for plotting data.
    %   simInfo - Simulation information structure, containing settings like the working directory.
    %   Optional:
    %     'savePlots' - Boolean to save the figure (default false)
    %     'showPlots' - Boolean to show the figure (default true)

    p = inputParser;
    addParameter(p, 'savePlots', true, @islogical);
    addParameter(p, 'showPlots', false, @islogical);
    parse(p, varargin{:});

    savePlots = p.Results.savePlots;
    showPlots = p.Results.showPlots;
    wdSim = simInfo.wdSim;
    T = simInfo.T;
    PSubs = sysInfo.PSubs_allT_vs_macroItr;
    PSubsCost = sysInfo.PSubsCost_allT_vs_macroItr;
    macroItr = 1:length(PSubs);  % Assuming the length of macro iterations is the same for both metrics

    % Create figure with controlled visibility
    f = figure('Visible', showPlots);

    % First y-axis (left side) for Substation Power
    % yyaxis left;
    % plot(macroItr, PSubs, '-o', 'Color', [0, 0.447, 0.741], 'MarkerFaceColor', [0, 0.447, 0.741], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Substation Power');
    % ylabel('Power [MW]', 'Color', [0, 0.447, 0.741]);
    % ylim([min(PSubs)*0.9, max(PSubs)*1.1]);

    % Second y-axis (right side) for Substation Cost
    % yyaxis right;
    plot(macroItr, PSubsCost, '-s', 'Color', [0.85, 0.325, 0.098], 'MarkerFaceColor', [0.85, 0.325, 0.098], 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Substation Power Cost');
    ylabel('Cost [\$]');
    ylim([min(PSubsCost)*0.9, max(PSubsCost)*1.1]);
    
    % hold on;

    % title('Substation Power and Cost Convergence');
    title('Convergence of Cost of Substation Power for Entire Horizon');
    xlabel('Macro Iteration Number');
    grid minor;
    legend('show', 'Location', 'northeast');

    % Adjust y-axis limits based on data
    % yyaxis left;
    % ylim([min(PSubs)*0.9, max(PSubs)*1.1]);
    % yyaxis right;
    % ylim([min(PSubsCost)*0.9, max(PSubsCost)*1.1]);

    % Save the figure if requested
    if savePlots
        folderPath = fullfile(wdSim, 'processedData', sysInfo.systemName, strcat('numAreas', num2str(sysInfo.numAreas)));
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        % filename = fullfile(folderPath, 'ObjectiveConvergenceCurves.png');
        ext = ".png";
        filename = fullfile(folderPath, strcat('ObjectiveConvergenceCurves_Horizon_', num2str(T), ext));
        saveas(f, filename);
    end

    hold off;
end
