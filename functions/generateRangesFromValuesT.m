function outputArrays = generateRangesFromValuesT(values, T)
    numVarsPerTimePeriod = sum(values);  % Total number of variables in one time step
    indexRangeBorders = cumsum([0, values]);
    numArrays = numel(indexRangeBorders) - 1;
    outputArrays = cell(1, numArrays);

    for i = 1:numArrays
        numIndices = indexRangeBorders(i+1) - indexRangeBorders(i);
        tempMatrix = zeros(T, numIndices);

        for t = 1:T
            offset = (t-1) * numVarsPerTimePeriod;
            tempMatrix(t, :) = (indexRangeBorders(i)+1 + offset):(indexRangeBorders(i+1) + offset);
        end

        outputArrays{i} = tempMatrix;
    end
end
