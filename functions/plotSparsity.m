function plotSparsity(Aeq, beq)
    % plotSparsity: Plot sparsity patterns of Aeq and beq with visually aligned rows.
    % 
    % Inputs:
    %   - Aeq: Sparse matrix
    %   - beq: Sparse vector

    % Define custom colors
    wineRed = [0.6, 0.2, 0.2];
    sexyGreen = [0, 0.6, 0.3];

    % Create the figure and set its size and properties
    figure('Color', 'white', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

    % Convert sparse matrices to full to capture the sparsity using scatter
    [iAeq, jAeq, vAeq] = find(Aeq);
    [ibeq, ~] = find(beq);

    % Determine positions for the axes objects
    axWidthAeq = 0.7;
    axWidthBeq = 0.2;
    axPosAeq = [0.05, 0.1, axWidthAeq, 0.8];
    axPosBeq = [0.78, 0.1, axWidthBeq, 0.8];

    % Plot Aeq using scatter. Positive values are wine red and negative ones are sexy green.
    ax1 = axes('Position', axPosAeq);
    scatter(jAeq(vAeq < 0), iAeq(vAeq < 0), 20, sexyGreen, 'filled');  % Negative values: sexy green
    hold on;
    scatter(jAeq(vAeq > 0), iAeq(vAeq > 0), 20, wineRed, 'filled');  % Positive values: wine red
    title('Sparsity pattern of Aeq', 'FontWeight', 'bold', 'FontSize', 14);
    axis equal tight;
    set(ax1, 'FontSize', 12, 'FontWeight', 'bold', 'YDir', 'reverse');
    hold off;

    % Plot beq using scatter with dark green color
    ax2 = axes('Position', axPosBeq);
    scatter(ones(size(ibeq)), ibeq, 20, [0, 0.5, 0], 'filled');  % dark green color
    title('Sparsity pattern of beq', 'FontWeight', 'bold', 'FontSize', 14);
    axis equal tight;
    set(ax2, 'FontSize', 12, 'FontWeight', 'bold', 'YDir', 'reverse');

    % Add a title for the entire figure using annotation
    annotation('textbox', [0, 0.93, 1, 0.07], ...
               'String', 'Sparsity Patterns of Aeq and beq', ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'FontSize', 16, ...
               'FontWeight', 'bold');

    drawnow;
end
