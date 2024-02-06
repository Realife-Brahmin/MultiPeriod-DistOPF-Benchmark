function sparseArray = sparseArrayFromDense(denseArray, totalElements, specifiedIndices)
    width = size(denseArray, 2);
    % Initialize an array of zeros for all elements
    sparseArray = zeros(totalElements, width);

    % Iterate over the specified indices
    for i = 1:length(specifiedIndices)
        % Assign the corresponding value from the dense array to the specified index
        index = specifiedIndices(i);
        sparseArray(index, 1:width) = denseArray(i, 1:width);
    end
end
