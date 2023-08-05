xSys = struct(); 
xSys.Names = ["Time Period"].';
xSys.VarNames = ["timePeriods"].';
xSys.Units = ["#"].';
xSys.FullNames = ["Discrete Time Period"].';
xSys.FigureNames = ["NA"].';
xSys.Legends = ["t"].';
xSys.yLabelNames = strcat(xSys.Legends, ' \thinspace [', xSys.Units, ']');
% Convert the xSys struct to a table
xSysTable = struct2table(xSys);

% Display the table
disp(xSysTable);
writetable(xSysTable, strcat("plottingSpecifications\", "xSysTable", ".csv"))