function outputArrays = generateRangesFromValuesT(values, T)
    numVarsPerTimePeriod = sum(values);  % Total number of variables in one time step
    indexRangeBorders = cumsum([0, values]);
    numArrays = numel(indexRangeBorders) - 1;
    outputArrays = cell(1, numArrays);
    
    for i = 1:numArrays
        tempRanges = cell(1, T);
        for t = 1:T
            offset = (t-1) * numVarsPerTimePeriod;
            tempRanges{t} = (indexRangeBorders(i)+1 + offset):(indexRangeBorders(i+1) + offset);
        end
        outputArrays{i} = tempRanges;
    end
end
