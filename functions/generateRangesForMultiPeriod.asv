% function eqnIndices = generateRangesForMultiPeriod(valuesPeriodsSet1, numPeriodsSet1, valuesPeriodsSet2, numPeriodsSet2)
%     % Calculate the total number of variables in each time step set
%     numVarsPerTimePeriod1 = sum(valuesPeriodsSet1);
%     numVarsPerTimePeriod2 = sum(valuesPeriodsSet2);
% 
%     % Calculate the cumulative sums to determine borders of index ranges
%     indexRangeBorders1 = cumsum([0, valuesPeriodsSet1]);
%     indexRangeBorders2 = cumsum([0, valuesPeriodsSet2]);
% 
%     numArrays1 = numel(indexRangeBorders1) - 1;
%     numArrays2 = numel(indexRangeBorders2) - 1;
% 
%     eqnIndices = cell(1, max(numArrays1, numArrays2));
% 
%     % Generate indices for normal periods (1 to numPeriodsSet1)
%     for i = 1:numArrays1
%         numIndices = indexRangeBorders1(i+1) - indexRangeBorders1(i);
%         tempMatrix = zeros(numIndices, numPeriodsSet1);
%         for t = 1:numPeriodsSet1
%             offset = (t-1) * numVarsPerTimePeriod1;
%             tempMatrix(:, t) = (indexRangeBorders1(i) + 1 + offset):(indexRangeBorders1(i+1) + offset);
%         end
%         eqnIndices{i} = tempMatrix;
%     end
% 
%     % Append indices for end periods (last numPeriodsSet2 periods)
%     for i = 1:numArrays2
%         numIndices = indexRangeBorders2(i+1) - indexRangeBorders2(i);
%         tempMatrix = zeros(numIndices, numPeriodsSet2);
%         for t = 1:numPeriodsSet2
%             offset = numPeriodsSet1 * numVarsPerTimePeriod1 + (t-1) * numVarsPerTimePeriod2;  % Correct offset for different vars count
%             tempMatrix(:, t) = (indexRangeBorders2(i) + 1 + offset):(indexRangeBorders2(i+1) + offset);
%         end
%         if i <= numel(eqnIndices)
%             eqnIndices{i} = [eqnIndices{i}, tempMatrix];  % Append to the existing matrix
%         else
%             eqnIndices{i} = tempMatrix;  % Handle case if more classes at end than in normal periods
%         end
%     end
% end
function eqnIndices = generateRangesForMultiPeriod(valuesPeriodsSet1, numPeriodsSet1, valuesPeriodsSet2, numPeriodsSet2)
    % Calculate the total number of variables for each period set
    numVarsPerTimePeriod1 = sum(valuesPeriodsSet1);
    numVarsPerTimePeriod2 = sum(valuesPeriodsSet2);

    % Calculate the cumulative sums to determine borders of index ranges for both sets
    indexRangeBorders1 = cumsum([0, valuesPeriodsSet1]);
    indexRangeBorders2 = cumsum([0, valuesPeriodsSet2]);

    % Determine the maximum number of arrays needed
    numArrays1 = numel(indexRangeBorders1) - 1;
    numArrays2 = numel(indexRangeBorders2) - 1;
    numArrays = max(numArrays1, numArrays2);

    eqnIndices = cell(1, numArrays);

    % Generate indices for Set1 (periods 1 to numPeriodsSet1)
    for i = 1:numArrays1
        tempIndices = [];
        for t = 1:numPeriodsSet1
            offset = (t-1) * numVarsPerTimePeriod1;
            tempIndices = [tempIndices, (indexRangeBorders1(i) + 1 + offset):(indexRangeBorders1(i+1) + offset)];
        end
        eqnIndices{i} = reshape(tempIndices, [], numPeriodsSet1);
    end

    % Append indices for Set2 (last numPeriodsSet2 periods)
    for i = 1:numArrays2
        tempIndices = [];
        for t = 1:numPeriodsSet2
            offset = numPeriodsSet1 * numVarsPerTimePeriod1 + (t-1) * numVarsPerTimePeriod2;
            tempIndices = [tempIndices, (indexRangeBorders2(i) + 1 + offset):(indexRangeBorders2(i+1) + offset)];
        end
        if i <= numArrays1 && size(eqnIndices{i}, 1) == length(tempIndices) / numPeriodsSet2
            eqnIndices{i} = [eqnIndices{i}, reshape(tempIndices, [], numPeriodsSet2)];
        else
            if iscell(eqnIndices{i})
                eqnIndices{i}{end+1} = reshape(tempIndices, [], numPeriodsSet2);
            else
                % Convert existing matrix to cell and append new data
                eqnIndices{i} = {eqnIndices{i}, reshape(tempIndices, [], numPeriodsSet2)};
            end
        end
    end
end

