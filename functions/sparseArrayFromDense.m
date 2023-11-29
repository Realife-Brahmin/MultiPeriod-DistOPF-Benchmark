function sparseArray = sparseArrayFromDense(denseArray, totalElements, specifiedIndices)
    % Initialize an array of zeros for all elements
    sparseArray = zeros(totalElements, 1);

    % Iterate over the specified indices
    for i = 1:length(specifiedIndices)
        % Assign the corresponding value from the dense array to the specified index
        index = specifiedIndices(i);
        sparseArray(index) = denseArray(i);
    end
end
