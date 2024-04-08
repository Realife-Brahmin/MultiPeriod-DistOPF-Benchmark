% Define a helper function for simultaneous console and file output
function print_and_write(fileID, locNum, formatSpec, varargin)
    fprintf(1, [num2str(locNum), '. ', formatSpec], varargin{:}); % Print to Command Window
    fprintf(fileID, [num2str(locNum), '. ', formatSpec], varargin{:}); % Write to file
end