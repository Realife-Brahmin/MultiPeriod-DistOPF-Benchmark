% function indices_v = excludeFirstElement(indices_vAll)
%     %EXCLUDEFIRSTELEMENT Exclude the first element of arrays in a cell array.
%     %
%     % Synopsis:
%     %   indices_v = excludeFirstElement(indices_vAll)
%     %
%     % Description:
%     %   Given a cell array of arrays, the function returns a cell array
%     %   where each array has its first element removed.
%     %
%     % Input:
%     %   - indices_vAll: Cell array of arrays.
%     %
%     % Output:
%     %   - indices_v: Modified cell array.
%     %
% 
%     indices_v = cellfun(@(x) x(2:end), indices_vAll, 'UniformOutput', false);
% end
function matrixReduced = excludeFirstElement(matrix2D)
    %EXCLUDEFIRSTELEMENT Exclude the first column of a 2D matrix.
    %
    % Synopsis:
    %   matrixReduced = excludeFirstElement(matrix2D)
    %
    % Description:
    %   Given a 2D matrix, the function returns another 2D matrix
    %   where the first column is removed.
    %
    % Input:
    %   - matrix2D: Input 2D matrix.
    %
    % Output:
    %   - matrixReduced: Modified 2D matrix.
    %

    matrixReduced = matrix2D(:, 2:end);
end


