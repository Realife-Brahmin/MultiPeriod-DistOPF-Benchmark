locNum = 0; % Initialize line number counter

% Starting block
disp('------------------------------------------------------------')
locNum = locNum + 1; disp([num2str(locNum) '. Machine ID: ', getenv("COMPUTERNAME")])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon Duration: ', num2str(T)])
locNum = locNum + 1; disp([num2str(locNum) '. Nature of Simulation: ', simNatureStringFull])
locNum = locNum + 1; disp([num2str(locNum) '. GED Configuruation: ', battstring])

for t = 1:T
    locNum = 0;
    disp('-----------------------------')
    disp(['Hour: ', num2str(t)])
    locNum = locNum + 1; disp([num2str(locNum) '. Line Loss: ', num2str(lineLoss_kW_1toT(t)),' kW'])
    locNum = locNum + 1; disp([num2str(locNum) '. Substation Power: ', num2str(substationPower_kW_1toT(t)),' kW + ', num2str(substationPower_kVAr_1toT(t)), ' kVAr'])
    locNum = locNum + 1; disp([num2str(locNum) '. Total Load: ', num2str(pLTotal_kW_1toT(t)), ' kW + ', num2str(qLTotal_kVAr_1toT(t)), ' kVAr'])
    locNum = locNum + 1; disp([num2str(locNum) '. Total Generation: ', num2str(pTotal_kW_1toT(t)), ' kW + ', num2str(qTotal_kVAr_1toT(t)), ' kVAr' ])
    locNum = locNum + 1; disp([num2str(locNum) '. Total PV Generation: ', num2str(pDTotal_kW_1toT(t)), ' kW + ', num2str(qDTotal_kVAr_1toT(t)), ' kVAr'])
    locNum = locNum + 1; disp([num2str(locNum) '. Total Battery Generation: ', num2str(PdcTotal_kW_1toT(t)), ' kW + ', num2str(qBTotal_kVAr_1toT(t)), ' kVAr'])
    locNum = locNum + 1; disp([num2str(locNum) '. Total Static Capacitor Reactive Power Generation: ', num2str(qCTotal_kVAr_1toT(t)), ' kVAr'])
    locNum = locNum + 1; disp([num2str(locNum) '. Substation Power Cost: $ ', num2str(genCost_dollars_1toT(t))])
    locNum = locNum + 1; disp([num2str(locNum) '. SCD Observed: ', num2str(P_scd_Total_kW_1toT(t)), ' kW'])
end

locNum = 0;

% Horizon summary block
disp('-----------------------------')
disp(['Hour: Full ', num2str(T), ' Hour Horizon'])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon Line Loss: ', num2str(lineLoss_kW_allT),' kW'])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon Total Substation Power: ', num2str(substationPower_kW_allT),' kW + ', num2str(substationPower_kVAr_allT), ' kVAr'])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon Total Load: ', num2str(pL_Total_kW_allT), ' kW + ', num2str(qL_Total_kVAr_allT), ' kVAr'])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon Total Generation: ', num2str(p_Total_kW_allT), ' kW + ', num2str(q_Total_kVAr_allT), ' kVAr' ])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon Total PV Generation: ', num2str(pD_Total_kW_allT), ' kW + ', num2str(qD_Total_kVAr_allT), ' kVAr'])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon Total Battery Generation: ', num2str(Pdc_Total_kW_allT), ' kW + ', num2str(qB_Total_kVAr_allT), ' kVAr'])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon Total Static Capacitor Reactive Power Generation: ', num2str(QC_Total_kVAr_allT), ' kVAr'])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon Total Substation Power Cost: $ ', num2str(genCost_dollars_allT)])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon Total SCD Observed: ', num2str(P_scd_Total_kW_allT), ' kW'])
locNum = locNum + 1; disp([num2str(locNum) '. Horizon-end Battery Energy Deviation from Reference: ', num2str(B_violation_abs_Total_kWh), ' kWh'])

locNum = 0;

% Simulation metadata block
disp('------------------------------------------------------------')
locNum = locNum + 1; disp([num2str(locNum) '. Number of Macro-Iterations: ', num2str(macroItr+1)])
locNum = locNum + 1; disp([num2str(locNum) '. Simulation Time: ', num2str(grandTotalTime), ' s'])
locNum = locNum + 1; disp([num2str(locNum) '. Time to solve with sequential (non-parallel) computation: ', num2str(time_if_serial), ' s'])
locNum = locNum + 1; disp([num2str(locNum) '. Time to solve if OPF computation parallelized: ', num2str(time_if_parallel), ' s'])
disp('------------------------------------------------------------')
