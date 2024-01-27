function outputArray = select_percentage_of_nz_elements(allBuses, DER_percent)
    % SELECT_PERCENTAGE_OF_NZ_ELEMENTS selects a percentage of non-zero elements
    % from a given array in a spatially even manner.
    %
    % Usage:
    %   outputArray = select_percentage_of_nz_elements(allBuses, DER_percent)
    %
    % Inputs:
    %   allBuses - A numeric array representing buses. Non-zero values indicate
    %              buses with DERs (Distributed Energy Resources).
    %   DER_percent - A scalar value representing the percentage of buses with DERs
    %                 to be selected. This is a value between 0 and 100.
    %
    % Outputs:
    %   outputArray - A numeric array of the same size as allBuses. Non-zero values
    %                 in outputArray indicate the selected buses based on the given
    %                 DER_percent.
    %
    % Description:
    %   The function selects a specified percentage of buses that have DERs (non-zero
    %   elements in allBuses). It does so by calculating evenly spaced indices using
    %   linspace, ensuring a spatially even distribution of selected buses. The
    %   outputArray contains the selected buses, with all other elements set to zero.
    %
    % Example:
    %   % Example usage with a hypothetical array of buses and a DER percent of 50
    %   allBuses = [1 0 2 0 3 4 0 5];
    %   DER_percent = 50;
    %   selectedBuses = select_percentage_of_nz_elements(allBuses, DER_percent);
    %
    % Notes:
    %   - The function uses the floor function on the calculated indices from linspace
    %     to ensure integer index values.
    %   - If DER_percent is set to 0, no buses will be selected, and the outputArray
    %     will be an array of zeros.
    %   - The ceil function is used to determine the number of elements to select,
    %     ensuring that at least one element is selected when DER_percent > 0.

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
