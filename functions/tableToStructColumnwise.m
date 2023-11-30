function s = tableToStructColumnwise(t)
    % Convert the table to a structure with each table column as a field
    varNames = t.Properties.VariableNames;
    for i = 1:length(varNames)
        s.(varNames{i}) = t.(varNames{i});
    end
end
