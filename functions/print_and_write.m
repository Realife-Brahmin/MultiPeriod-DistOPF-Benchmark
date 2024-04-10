% Define a helper function for simultaneous console and file output
function print_and_write(fileID, locNum, formatSpec, varargin)
    locNumStr = sprintf('%02d', locNum);  % Use sprintf to format the number with leading zeros
    fprintf(1, [locNumStr, '. ', formatSpec], varargin{:}); % Print to Command Window
    fprintf(fileID, [locNumStr, '. ', formatSpec], varargin{:}); % Write to file
end
