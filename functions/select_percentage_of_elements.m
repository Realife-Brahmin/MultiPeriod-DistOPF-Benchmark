function selected_elements = select_percentage_of_elements(array, percentage)
    % Ensure the percentage is between 0 and 100
    if percentage < 0 || percentage > 100
        error('Percentage must be between 0 and 100');
    end
    
    % Calculate the number of elements to select
    num_elements = numel(array);
    num_to_select = ceil(num_elements * (percentage / 100));
    
    % Calculate the starting index to extract the middle percentage
    startIndex = ceil((num_elements - num_to_select + 1) / 2);
    
    % Extract the middle percentage of elements
    endIndex = startIndex + num_to_select - 1;  % Ensure that endIndex does not exceed the array length
    endIndex = min(endIndex, num_elements);     % Correct endIndex if necessary
    
    selected_elements = array(startIndex:endIndex);
end
