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

    % Plot Aeq using scatter
    ax1 = axes('Position', axPosAeq);
    scatter(jAeq(vAeq < 0), iAeq(vAeq < 0), 20, sexyGreen, 'filled');
    hold on;
    scatter(jAeq(vAeq > 0), iAeq(vAeq > 0), 20, wineRed, 'filled');
    title('Sparsity pattern of Aeq', 'FontWeight', 'bold', 'FontSize', 14);
    axis equal tight;
    set(ax1, 'FontSize', 12, 'FontWeight', 'bold', 'YDir', 'reverse');
    
    % Set major and minor grid properties for Aeq
    ax1.GridAlpha = 0.3;
    ax1.MinorGridAlpha = 0.2;
    ax1.XAxis.MinorTickValues = 0.5:1:(size(Aeq,2) + 0.5);
    ax1.YAxis.MinorTickValues = 0.5:1:(size(Aeq,1) + 0.5);
    grid(ax1, 'on');
    grid(ax1, 'minor');
    hold off;

    % Plot beq using scatter
    ax2 = axes('Position', axPosBeq);
    scatter(ones(size(ibeq)), ibeq, 20, [0, 0.5, 0], 'filled');
    title('Sparsity pattern of beq', 'FontWeight', 'bold', 'FontSize', 14);
    axis equal tight;
    set(ax2, 'FontSize', 12, 'FontWeight', 'bold', 'YDir', 'reverse');
    
    % Set major and minor grid properties for beq
    ax2.GridAlpha = 0.3;
    ax2.MinorGridAlpha = 0.2;
    ax2.YAxis.MinorTickValues = 0.5:1:(size(beq,1) + 0.5);
    grid(ax2, 'on');
    grid(ax2, 'minor');

    % Add a title for the entire figure using annotation
    annotation('textbox', [0, 0.93, 1, 0.07], ...
               'String', 'Sparsity Patterns of Aeq and beq', ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'FontSize', 16, ...
               'FontWeight', 'bold');

    drawnow;
end
