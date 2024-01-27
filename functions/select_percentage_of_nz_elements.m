function outputArray = select_percentage_of_nz_elements(allBuses, DER_percent)
    % Extract buses with DERs (non-zero elements)
    busesWithDERs = find(allBuses);

    % Calculate the number of elements to select
    numToSelect = ceil(length(busesWithDERs) * DER_percent / 100);

    % Use linspace to generate evenly spaced indices
    indices = floor(linspace(1, length(busesWithDERs), numToSelect));

    % Create the output array, initializing with zeros
    outputArray = zeros(size(allBuses));

    % If there are buses to select
    if numToSelect > 0
        % Select buses based on the indices from busesWithDERs
        selectedBuses = busesWithDERs(indices);
        outputArray(selectedBuses) = allBuses(selectedBuses);
    end
end
