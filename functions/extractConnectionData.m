function [CBMatrix, CBTable] = extractConnectionData(systemDataFolder, varargin)
    % EXTRACTCONNECTIONDATA - Extract connection data from a system data folder
    %
    %   extractConnectionData(systemDataFolder) extracts connection data
    %   from the specified system data folder.
    %
    %   extractConnectionData(systemDataFolder, 'verbose', true) extracts
    %   connection data with verbose output enabled.
    %
    %   Input arguments:
    %   - systemDataFolder: The path to the system data folder.
    %
    %   Optional arguments:
    %   - 'verbose': A logical value indicating whether table should
    %     be displayed (default: false).
    
    % Default values
    verbose = false;
    
    % Parse optional arguments
    for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'verbose')
            verbose = varargin{i+1};
        end
    end

    connectionBusesFilename = strcat(systemDataFolder, "connectionBuses.txt");
    CBTable0 = readtable(connectionBusesFilename);
    CBMatrix = table2array(CBTable0);
    
    connectionBuses_FullFilename = strcat(systemDataFolder, "connectionBuses_Full.txt");
    CBTable = readtable(connectionBuses_FullFilename);
    CBTable.Properties.VariableNames(3:6) = ["conBus_parentAreaFrom", "conBus_parentAreaTo", "conBus_childAreaFrom", "conBus_childAreaTo"];
    CBTable.S_childArea = zeros(size(CBTable, 1), 1);
    CBTable.v2_parentArea = 1.03^2*ones(size(CBTable, 1), 1);
    mydisplay(verbose, CBTable);
 
end