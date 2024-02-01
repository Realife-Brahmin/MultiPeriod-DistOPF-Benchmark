function sysInfoTruncated = truncateSysInfo(sysInfo, macroItr)
    % Truncates specified fields of the sysInfo structure to include only non-zero instances
    % up to macroItr+1, while keeping other fields unchanged.

    % List of fields to truncate
    fieldsToTruncate = {'PLoss_allT_vs_macroItr', 'PLoss_1toT_vs_macroItr', ...
                        'PSubs_allT_vs_macroItr', 'PSubs_1toT_vs_macroItr', ...
                        'PSubsCost_allT_vs_macroItr', 'PSubsCost_1toT_vs_macroItr'};

    % Copy the entire sysInfo structure
    sysInfoTruncated = sysInfo;

    for i = 1:length(fieldsToTruncate)
        field = fieldsToTruncate{i};

        % Check if the field exists in sysInfo
        if isfield(sysInfo, field)
            data = sysInfo.(field);

            % Determine the size of the data
            [rows, cols] = size(data);

            % Truncate the data
            if cols == 1  % For column vectors
                truncatedData = data(1:min(rows, macroItr + 1));
            else  % For matrices
                truncatedData = data(:, 1:min(cols, macroItr + 1));
            end

            % Assign the truncated data to the new structure
            sysInfoTruncated.(field) = truncatedData;
        end
    end
end
