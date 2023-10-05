% plottingSpecificationFolder = strcat("plottingSpecifications", filesep);

% tableNames = {'xSys', 'ySys', 'yArs', 'yNodes'};
% for i = 1:length(tableNames)
%     tableType = tableNames{i};
%     tableName = strcat(tableType, "Table");
%     table = readtable(fullfile(plottingSpecificationFolder, strcat(tableName, ".csv")));
%     eval([tableType, ' = table2struct(table, "ToScalar", true);']);
% end

% timePeriods = transpose(1:T);