% Assuming sysInfo, numAreas, T, and noBatteries are already defined in your workspace

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

% Initialize line number counter for the text file
fileLocNum = 0;

% Starting block
fprintf(fileID, '------------------------------------------------------------\n');
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Machine ID: %s\n', fileLocNum, getenv("COMPUTERNAME"));
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon Duration: %d\n', fileLocNum, T);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Nature of Simulation: %s\n', fileLocNum, simNatureStringFull);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. GED Configuration: %s\n', fileLocNum, battstring);

% Write the results to the file for each time step
for t = 1:T
    fprintf(fileID, '-----------------------------\n');
    fprintf(fileID, 'Hour: %d\n', t);
    fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Line Loss: %.2f kW\n', fileLocNum, lineLoss_kW_1toT(t));
    fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Substation Power: %.2f kW + %.2f kVAr\n', fileLocNum, substationPower_kW_1toT(t), substationPower_kVAr_1toT(t));
    fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Total Load: %.2f kW + %.2f kVAr\n', fileLocNum, pLTotal_kW_1toT(t), qLTotal_kVAr_1toT(t));
    fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Total Generation: %.2f kW + %.2f kVAr\n', fileLocNum, pTotal_kW_1toT(t), qTotal_kVAr_1toT(t));
    fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Total PV Generation: %.2f kW + %.2f kVAr\n', fileLocNum, pDTotal_kW_1toT(t), qDTotal_kVAr_1toT(t));
    fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Total Battery Generation: %.2f kW + %.2f kVAr\n', fileLocNum, PdcTotal_kW_1toT(t), qBTotal_kVAr_1toT(t));
    fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Total Static Capacitor Reactive Power Generation: %.2f kVAr\n', fileLocNum, qCTotal_kVAr_1toT(t));
    fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Substation Power Cost: $%.3f\n', fileLocNum, genCost_dollars_1toT(t));
    fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. SCD Observed: %.2f kW\n', fileLocNum, P_scd_Total_kW_1toT(t));
end

% Horizon summary block
fprintf(fileID, '-----------------------------\n');
fprintf(fileID, 'Hour: Full %d Hour Horizon\n', T);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon Line Loss: %.2f kW\n', fileLocNum, lineLoss_kW_allT);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon Total Substation Power: %.2f kW + %.2f kVAr\n', fileLocNum, substationPower_kW_allT, substationPower_kVAr_allT);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon Total Load: %.2f kW + %.2f kVAr\n', fileLocNum, pL_Total_kW_allT, qL_Total_kVAr_allT);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon Total Generation: %.2f kW + %.2f kVAr\n', fileLocNum, p_Total_kW_allT, q_Total_kVAr_allT);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon Total PV Generation: %.2f kW + %.2f kVAr\n', fileLocNum, pD_Total_kW_allT, qD_Total_kVAr_allT);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon Total Battery Generation: %.2f kW + %.2f kVAr\n', fileLocNum, Pdc_Total_kW_allT, qB_Total_kVAr_allT);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon Total Static Capacitor Reactive Power Generation: %.2f kVAr\n', fileLocNum, QC_Total_kVAr_allT);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon Total Substation Power Cost: $%.3f\n', fileLocNum, genCost_dollars_allT);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon Total SCD Observed: %.2f kW\n', fileLocNum, P_scd_Total_kW_allT);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Horizon-end Battery Energy Deviation from Reference: %.2f kWh\n', fileLocNum, B_violation_abs_Total_kWh);

% Simulation metadata block
fprintf(fileID, '------------------------------------------------------------\n');
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Number of Macro-Iterations: %d\n', fileLocNum, macroItr+1);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Simulation Time: %.2f s\n', fileLocNum, grandTotalTime);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Time to solve with sequential (non-parallel) computation: %.2f s\n', fileLocNum, time_if_serial);
fileLocNum = fileLocNum + 1; fprintf(fileID, '%d. Time to solve if OPF computation parallelized: %.2f s\n', fileLocNum, time_if_parallel);
fprintf(fileID, '------------------------------------------------------------\n');

% Close the file
fclose(fileID);

% Display a message indicating that the results were saved
disp(['Results saved to ', systemsSolutionName_fval]);
