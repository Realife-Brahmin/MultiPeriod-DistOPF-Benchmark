ySys = struct(); 
ySys.Names = ["PLoss", "QLoss", "PSubs", "QSubs", "TimeToSolve"].';
ySys.VarNames = ["PLosses", "QLosses", "PSubs", "QSubs", "solutionTimes"].';
ySys.Units = ["kW", "kVAr", "kW", "kVAr", "s"].';
ySys.FullNames = ["Line Losses", "Line Reactive Losses", "Substation Power", "Substation Reactive Power", "Solution Times"].';
ySys.FigureNames = ["lossesFigure", "reactiveLossesFigure", "substationPowerFigure", "substationReactiveFigure", "solutionTimesFigure"].';
ySys.legends = ["P_{Loss}", "Q_{Loss}", "P_{Subs}", "Q_{Subs}", "t"].';
ySys.yLabelNames = strcat(ySys.legends, ' \thinspace [', ySys.Units, ']');
% Convert the ySys struct to a table
ySysTable = struct2table(ySys);

% Display the table
disp(ySysTable);

writetable(ySysTable, strcat("plottingSpecifications\", "ySysTable", ".csv"))