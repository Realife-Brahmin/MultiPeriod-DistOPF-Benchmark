function middleData = choose_middle_T_points(data, T)
    % Your data column vector
    % data = rand(24, 1); % Example data
    
    % The number of hours for your simulation
    % T = 10; % Example value
    
    % Calculate the starting index to extract the middle T points
    startIndex = ceil((length(data) - T + 1) / 2);
    
    % Extract the middle T points
    middleData = data(startIndex:startIndex+T-1);
end