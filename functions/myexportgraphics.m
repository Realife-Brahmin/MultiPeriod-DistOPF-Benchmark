function myexportgraphics(saveFlag, varargin)
% MYEXPORTGRAPHICS Exports a figure using exportgraphics if saveFlag is true.
%   If saveFlag is true, this function uses exportgraphics to export the
%   current figure in the desired format specified by the additional
%   arguments. If saveFlag is false, the figure is displayed on the screen
%   without exporting.

% Check if saveFlag is true
if saveFlag
    % Use exportgraphics to save the figure
    exportgraphics(varargin{:});
else
    % Display the figure on the screen without saving
    disp('Figure not saved.');
end

end