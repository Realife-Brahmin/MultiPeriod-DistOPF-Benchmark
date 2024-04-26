function costArray = generateCostProfileBasedOnLoads_OnPeakOffPeak(lambdaVals, lo, hi, peakHoursFraction)
    % Calculate the number of peak hours, ensuring at least one peak hour is counted
    T = length(lambdaVals);
    num_peak_hours = max(1, floor(peakHoursFraction * T));

    % Create the output array initialized to the low value
    costArray = lo * ones(T, 1);

    % Find indices for the highest values using sorting
    [~, sortedIndices] = sort(lambdaVals, 'descend');
    
    % Assign the high value to the peak hour indices
    costArray(sortedIndices(1:num_peak_hours)) = hi;
end
