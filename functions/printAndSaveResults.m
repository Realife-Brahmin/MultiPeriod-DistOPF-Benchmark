% Assuming sysInfo, numAreas, T, and other necessary variables are already defined in your workspace

% Construct the folder and file names
folderName = fullfile('processedData', sysInfo.systemName, strcat('numAreas_', num2str(numAreas), "/"));
prefixName = strcat(folderName, strcat('Horizon_', num2str(T)));

% Construct the full path for the TXT file
systemsSolutionName_fval = strcat(prefixName, "_", getenv('COMPUTERNAME'), '_results_', battstring, '.txt');

% Open the file for writing
fileID = fopen(systemsSolutionName_fval, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Failed to open file for writing.');
end

% Initialize line number counter
locNum = 0;

% Starting block
fprintf(fileID, '------------------------------------------------------------\n');
fprintf(1, '------------------------------------------------------------\n');
locNum = locNum + 1; print_and_write(fileID, locNum, 'Machine ID: %s\n', getenv("COMPUTERNAME"));
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Duration: %d\n', T);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Nature of Simulation: %s\n', simNatureStringFull);
locNum = locNum + 1; print_and_write(fileID, locNum, 'GED Configuration: %s\n', battstring);

% Loop for hourly details
% Print and write loop for hourly details
for t = 1:T
    locNum = 0;
    fprintf(fileID, '-----------------------------\n');
    fprintf(1, '-----------------------------\n');
    fprintf(fileID, 'Hour: %d\n', t);
    fprintf(1, 'Hour: %d\n', t);
    locNum = locNum + 1; print_and_write(fileID, locNum, 'Line Loss: %.2f kW\n', lineLoss_kW_1toT(t));
    locNum = locNum + 1; print_and_write(fileID, locNum, 'Substation Power: %.2f kW + %.2f kVAr\n', substationPower_kW_1toT(t), substationPower_kVAr_1toT(t));
    locNum = locNum + 1; print_and_write(fileID, locNum, 'Total Load: %.2f kW + %.2f kVAr\n', pLTotal_kW_1toT(t), qLTotal_kVAr_1toT(t));
    locNum = locNum + 1; print_and_write(fileID, locNum, 'Total Generation: %.2f kW + %.2f kVAr\n', pTotal_kW_1toT(t), qTotal_kVAr_1toT(t));
    locNum = locNum + 1; print_and_write(fileID, locNum, 'Total PV Generation: %.2f kW + %.2f kVAr\n', pDTotal_kW_1toT(t), qDTotal_kVAr_1toT(t));
    locNum = locNum + 1; print_and_write(fileID, locNum, 'Total PV Generation: %.2f kW + %.2f kVAr\n', pDTotal_kW_1toT(t), qDTotal_kVAr_1toT(t));
    locNum = locNum + 1; print_and_write(fileID, locNum, 'Total Battery Generation: %.2f kW + %.2f kVAr\n', PdcTotal_kW_1toT(t), qBTotal_kVAr_1toT(t));
    locNum = locNum + 1; print_and_write(fileID, locNum, 'Total Battery Flexibility Utilized: %.2f kW + %.2f kVAr\n', Pbatt_abs_Total_kW_1toT(t), qB_abs_Total_kVAr_1toT(t));
    locNum = locNum + 1; print_and_write(fileID, locNum, 'Total Static Capacitor Reactive Power Generation: %.2f kVAr\n', qCTotal_kVAr_1toT(t));
    locNum = locNum + 1; print_and_write(fileID, locNum, 'Substation Power Cost: $%.2f\n', genCost_dollars_1toT(t));
    locNum = locNum + 1; print_and_write(fileID, locNum, 'SCD Observed: %.2f kW\n', P_scd_Total_kW_1toT(t));
    % Add more metrics as necessary...
end

locNum = 0;
% Horizon summary block
fprintf(fileID, '-----------------------------\n');
fprintf(1, '-----------------------------\n');
fprintf(fileID, 'Hour: Full %d Hour Horizon\n', T);
fprintf(1, 'Hour: Full %d Hour Horizon\n', T);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Line Loss: %.2f kW\n', lineLoss_kW_allT);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Total Substation Power: %.2f kW + %.2f kVAr\n', substationPower_kW_allT, substationPower_kVAr_allT);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Total Load: %.2f kW + %.2f kVAr\n', pL_Total_kW_allT, qL_Total_kVAr_allT);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Total Generation: %.2f kW + %.2f kVAr\n', p_Total_kW_allT, q_Total_kVAr_allT);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Total PV Generation: %.2f kW + %.2f kVAr\n', pD_Total_kW_allT, qD_Total_kVAr_allT);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Total Battery Generation: %.2f kW + %.2f kVAr\n', Pdc_Total_kW_allT, qB_Total_kVAr_allT);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Total Battery Flexibility Utilized: %.2f kW + %.2f kVAr\n', Pbatt_abs_Total_kW_allT, qB_abs_Total_kVAr_allT);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Total Static Capacitor Reactive Power Generation: %.2f kVAr\n', QC_Total_kVAr_allT);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Total Substation Power Cost: $%.2f\n', genCost_dollars_allT);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon Total SCD Observed: %.2f kW\n', P_scd_Total_kW_allT);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Horizon-end Battery Energy Deviation from Reference: %.2f kWh\n', B_violation_abs_Total_kWh);

% locNum = 0;
% Simulation metadata block
fprintf(fileID, '------------------------------------------------------------\n');
fprintf(1, '------------------------------------------------------------\n');
locNum = locNum + 1; print_and_write(fileID, locNum, 'Number of Macro-Iterations: %d\n', macroItr+1);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Simulation Time: %.2f s\n', grandTotalTime);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Time to solve with sequential (non-parallel) computation: %.2f s\n', time_if_serial);
locNum = locNum + 1; print_and_write(fileID, locNum, 'Time to solve if OPF computation parallelized: %.2f s\n', time_if_parallel);
fprintf(fileID, '------------------------------------------------------------\n');
fprintf(1, '------------------------------------------------------------\n');

% Close the file
fclose(fileID);

% Display a message indicating that the results were saved
disp(strcat("Results saved to ", systemsSolutionName_fval));
