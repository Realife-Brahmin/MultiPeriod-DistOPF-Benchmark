% function result = getIndicesT(cellArray, index)
%     %extractFromCells Extracts the i-th element from each cell of the input cell array.
%     %
%     % Usage:
%     %   result = getIndicesT(cellArray, index)
%     %
%     % Input:
%     %   cellArray: A 1xN cell array where each cell contains an array of numbers.
%     %   index: The position of the element to extract from each array inside the cells.
%     %
%     % Output:
%     %   result: A 1xN array containing the extracted elements.
%     %
% 
%     result = cellfun(@(x) x(index), cellArray);
% end
function result = getIndicesT(matrix2D, col)
    %GETINDICEST Extracts the col-th element from each row of the input 2D array.
    %
    % Usage:
    %   result = getIndicesT(matrix2D, col)
    %
    % Input:
    %   matrix2D: A TxN array.
    %   col: The column number to extract from each row.
    %
    % Output:
    %   result: A Tx1 array containing the extracted elements.
    %

    result = matrix2D(:, col);
end

