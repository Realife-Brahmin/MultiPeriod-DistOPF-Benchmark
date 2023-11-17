function outputArray = select_percentage_of_nz_elements(inputArray, DER_percent)
    % Find the indices of non-zero elements
    nzIndices = find(inputArray ~= 0);
    % Calculate the number of non-zero elements to keep
    nzToKeep = ceil(length(nzIndices) * DER_percent / 100);
    % If nzToKeep is zero, we don't want to keep any elements
    if nzToKeep < 1
        outputArray = zeros(size(inputArray)); % Return an array of zeros
        return;
    end
    
    % Randomly select the indices of non-zero elements to keep
    selectedIndices = randsample(nzIndices, nzToKeep);
    
    % Create the output array, initializing with zeros
    outputArray = zeros(size(inputArray));
    % Place the selected non-zero elements in the output array
    outputArray(selectedIndices) = inputArray(selectedIndices);
end
