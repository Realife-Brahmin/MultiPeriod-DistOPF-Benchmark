function moveBatteryMonitorFiles(circuitName, busesWithBatts_Area, nBatt_Area, batteryMonitorFolder, varargin)
    % moveBatteryMonitorFiles - Move battery monitor CSV files to a specified folder
    %
    % Syntax:  
    %     moveBatteryMonitorFiles(circuitName, busesWithBatts_Area, nBatt_Area, batteryMonitorFolder, varargin)
    %
    % Description:
    %     moveBatteryMonitorFiles moves the specified battery monitor CSV files from the current directory 
    %     to a specified folder. It uses the bus numbers from the OpenDSS simulation associated with batteries 
    %     to identify the correct CSV files. Each file is named according to the convention: 
    %     <circuitName>_Mon_<monitorName>_1.csv, where <monitorName> is derived from the bus number.
    %
    % Inputs:
    %     circuitName          - String. Name of the OpenDSS circuit.
    %     busesWithBatts_Area  - Array. Bus numbers that have associated battery monitors.
    %     nBatt_Area           - Integer. Number of buses with batteries.
    %     batteryMonitorFolder - String. Path to the folder where CSV files will be moved.
    %     varargin             - (Optional) Name-Value pairs.
    %                           'verbose' - Logical. Set to true for display of operation details. Default is false.
    %
    % Outputs:
    %     None. The function moves files and optionally displays information about the operations.
    %
    % Examples:
    %     % Example 1: Move monitor files with verbose output
    %     moveBatteryMonitorFiles('MyCircuit', [1, 2, 3], 3, 'path/to/folder', 'verbose', true);
    %
    %     % Example 2: Move monitor files without verbose output
    %     moveBatteryMonitorFiles('MyCircuit', [1, 2, 3], 3, 'path/to/folder');
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: OTHER_FUNCTION_NAME
    %

    % Create an input parser
    p = inputParser;
    
    % Define default value for verbose
    defaultVerbose = false;
    
    % Add the verbose parameter to the parser
    addParameter(p, 'verbose', defaultVerbose, @islogical);
    
    % Parse the input arguments
    parse(p, varargin{:});
    
    % Retrieve the verbose value
    verbose = p.Results.verbose;

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
            
            % Display information if verbose is true
            if verbose
                disp(['Moved: ', csvFileName, ' to ', newFilePath]);
            end
        else
            if verbose
                disp(['File not found: ', csvFileName]);
            end
        end
    end
end
