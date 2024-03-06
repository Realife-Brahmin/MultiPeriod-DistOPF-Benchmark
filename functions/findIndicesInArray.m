function indices = findIndicesInArray(fullArray, subArray)
    
    % Preallocate the indices vector as a column vector with zeros, size equal to length of subArray
    indices = zeros(length(subArray), 1);
    
    % Loop through each element of the subArray
    for i = 1:length(subArray)
        % Find the index of the current element of subArray in fullArray
        idx = find(fullArray == subArray(i));
        
        % Assign the found index to the corresponding position in the indices vector
        indices(i) = idx;
    end
end
