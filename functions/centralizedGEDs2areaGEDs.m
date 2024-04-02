function P_GED_Area0 = centralizedGEDs2areaGEDs(busDataTable_Area, busesWithGEDs_Area0)
    % Initialize P_GED_Area0 with zeros based on the number of buses in the area
    P_GED_Area0 = zeros(size(busDataTable_Area.bus));
    
    % Iterate over each bus in the area
    for idx = 1:length(busDataTable_Area.bus)
        % Check if the current bus is in the list of buses with DERs and active
        if ismember(busDataTable_Area.bus(idx), busesWithGEDs_Area0)
            % Assign the P_der value from the busDataTable_Area if the bus is active
            % Note that by using `busDataTable_Area.P_der`, I'm still using
            % the original DER bus locations as my reference for all my GED
            % assignments
            P_GED_Area0(idx) = busDataTable_Area.P_der(idx);
        end
        % Note: If the bus is not in busesWithGEDs_Area0 or not active, P_GED_Area0(idx) remains zero
    end
end
