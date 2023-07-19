function outputArrays = generateRangesFromValues(values)
    indexRangeBorders = cumsum([0, values]);
    numArrays = numel(indexRangeBorders) - 1;
    outputArrays = cell(1, numArrays);
    
    for i = 1:numArrays
        outputArrays{i} = indexRangeBorders(i)+1:indexRangeBorders(i+1);
    end
end
