function result = getIndicesT(cellArray, index)
    %extractFromCells Extracts the i-th element from each cell of the input cell array.
    %
    % Usage:
    %   result = extractFromCells(cellArray, index)
    %
    % Input:
    %   cellArray: A 1xN cell array where each cell contains an array of numbers.
    %   index: The position of the element to extract from each array inside the cells.
    %
    % Output:
    %   result: A 1xN array containing the extracted elements.
    %

    result = cellfun(@(x) x(index), cellArray);
end
