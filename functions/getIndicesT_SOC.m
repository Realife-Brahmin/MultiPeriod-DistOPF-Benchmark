function result = getIndicesT_SOC(indicesSOC, batteryIndex)
    %GETINDICEST_SOC Extracts indices across periods with variable equation counts.
    %
    % Usage:
    %   result = getIndicesT_SOC(indicesSOC, batteryIndex)
    %
    % Input:
    %   indicesSOC: A cell array, where the first cell is a matrix for periods 1 to T-1,
    %               and the second cell is a vector for period T with double the entries per battery.
    %   batteryIndex: The index of the battery for which indices are required.
    %
    % Output:
    %   result: An array containing all indices for the specified battery across all periods.

    result = [];  % Initialize the result array

    % Handle periods 1 to T-1
    if ismatrix(indicesSOC{1})
        result = [result, indicesSOC{1}(batteryIndex, :)];
    end

    % Handle period T, assuming double the entries per battery
    if isvector(indicesSOC{2})
        numBatteries = length(indicesSOC{2}) / 2;  % There are double the entries, so divide by 4
        entriesPerBattery = length(indicesSOC{2}) / numBatteries;  % Calculate entries per battery
        baseIndex = (batteryIndex - 1) * entriesPerBattery + 1;
        % keyboard;
        result = [result, transpose(indicesSOC{2}(baseIndex:baseIndex + entriesPerBattery - 1))];
    end

    result = result';  % Transpose to match expected output format
end
