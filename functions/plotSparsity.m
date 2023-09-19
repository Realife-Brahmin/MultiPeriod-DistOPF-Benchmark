function plotSparsity(Aeq, beq)
    % plotSparsity: Plot sparsity patterns of Aeq and beq with visually aligned rows.
    % 
    % Inputs:
    %   - Aeq: Sparse matrix
    %   - beq: Sparse vector

    % Create the figure and set its size and properties
    figure('Color', 'white', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.5]);

    % Convert sparse matrices to full to capture the sparsity using scatter
    [iAeq, jAeq] = find(Aeq);
    [ibeq, ~] = find(beq);

    % Determine positions for the axes objects
    axWidthAeq = 0.6;
    axWidthBeq = 0.2;
    axHeight = size(Aeq, 1) / (size(Aeq, 1) + size(beq, 1));
    axPosAeq = [0.1, 0.1, axWidthAeq, axHeight];
    axPosBeq = [0.8, 0.1, axWidthBeq, axHeight];

    % Plot Aeq using scatter with wine red color and adjusted marker size
    ax1 = axes('Position', axPosAeq);
    scatter(jAeq, -iAeq, 20, [0.5, 0, 0.5], 'filled');  % wine red color
    title('Sparsity pattern of Aeq', 'FontWeight', 'bold', 'FontSize', 14);
    axis equal tight;
    grid on;
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');

    % Plot beq using scatter with dark green color
    ax2 = axes('Position', axPosBeq);
    scatter(ones(size(ibeq)), -ibeq, 20, [0, 0.5, 0], 'filled');  % dark green color
    title('Sparsity pattern of beq', 'FontWeight', 'bold', 'FontSize', 14);
    axis equal tight;
    grid on;
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');

    % Add a title for the entire figure using annotation
    annotation('textbox', [0, 0.9, 1, 0.1], ...
               'String', 'Sparsity Patterns of Aeq and beq', ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'FontSize', 16, ...
               'FontWeight', 'bold');

    drawnow;
end
