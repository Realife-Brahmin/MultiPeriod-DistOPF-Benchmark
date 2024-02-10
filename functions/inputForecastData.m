function [pvCoeffVals, lambdaVals, S_to_P_ratio_PV, S_to_P_ratio_Batt, costArray] = inputForecastData(rawDataFolder, T, globalPVCoeff)
    % Input: 
    % rawDataFolder - Folder containing the CSV files
    % T - Parameter for the choose_middle_T_points function
    % globalPVCoeff - Global PV Coefficient

    % Define the file extension
    ext = ".csv";

    % Construct file paths
    LoadShapePV24_addr = strcat(rawDataFolder, filesep, "LoadShapePV24", ext);
    LoadShape24_addr = strcat(rawDataFolder, filesep, "LoadShape24", ext);

    % Read the tables
    LoadShapePV24_table = readtable(LoadShapePV24_addr);
    LoadShape24_table = readtable(LoadShape24_addr);

    % Convert tables to arrays, assuming values are in the second column
    LoadShape24 = table2array(LoadShape24_table(:, 2));
    LoadShapePV24 = table2array(LoadShapePV24_table(:, 2));

    % Apply global PV Coefficient and choose middle T points
    pvCoeffVals = globalPVCoeff * choose_middle_T_points(LoadShapePV24, T);
    lambdaVals = choose_middle_T_points(LoadShape24, T);

    % Calculate S_to_P_ratio for PV and Battery
    S_to_P_ratio_PV = 1.2 * (min(1, max(pvCoeffVals)));
    S_to_P_ratio_Batt = S_to_P_ratio_PV;

    % Sample cost array and choosing middle T points
    costArray0 = [1.6, 3.8, 1.7, 1.8, 2.1, 2.2, 2.4, 2.9, 4.4, 5.1, 4.3, 3.5, 2.8, 2.8, 2.3, 2.0, 2.3, 4.1, 3.4, 6.6, 4.7, 4.7, 5.0, 2.6]';
    costArray = choose_middle_T_points(costArray0, T);

end