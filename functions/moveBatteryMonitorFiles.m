function moveBatteryMonitorFiles(circuitName, busesWithBatts_Area, nBatt_Area, batteryMonitorFolder)
    % Iterating through each bus with a battery
    for i = 1:nBatt_Area
        busNum = busesWithBatts_Area(i);
        monitorName = ['Battery', num2str(busNum), '_states']; % Construct the monitor name
        csvFileName = [circuitName, '_Mon_', monitorName, '_1.csv']; % Construct the CSV file name

        % Check if the file exists in the current directory
        if exist(csvFileName, 'file')
            % Construct the new file path
            newFilePath = fullfile(batteryMonitorFolder, csvFileName);

            % Move the file to the new location
            movefile(csvFileName, newFilePath);
            disp(['Moved: ', csvFileName, ' to ', newFilePath]);
        else
            disp(['File not found: ', csvFileName]);
        end
    end
end
