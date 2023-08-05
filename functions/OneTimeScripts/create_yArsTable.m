yArs = struct();
yArs.Names = ["P12", "Q12", "V1", "PSave"].';
yArs.Units = ["kW", "kVAr", "pu", "%"].';
yArs.VarNames = ["P12s", "Q12s", "V1s", "PSaves"].';
yArs.FullNames = ["Boundary Real Powers", "Boundary Reactive Powers", "Boundary Voltages", "Saved Power Percentage"].';
yArs.FigureNames = ["boundaryRealPowersFigure", "boundaryReactivePowersFigure", "boundaryVoltagesFigure", "savedPowerPercentage"].';
yArs.Legends = ["V_1", "P_{12}", "Q_{12}", "P_{Save}"].';
yArs.yLabelNames = strcat(yArs.Legends, " \thinspace [", yArs.Units, "]");
yArsTable = struct2table(yArs);

% Display the table
disp(yArsTable);
writetable(yArsTable, strcat("plottingSpecifications\", "yArsTable", ".csv"))