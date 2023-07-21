function postProcessResults(N, numAreas, N_Areas, isRoot, v2, actualBusNums, microIterationLosses, itr, time_dist, kVA_B, displaySimulationResultPlots, R_max, saveSimulationResultPlots, systemName, S_child, start, saveSimulationResults)
%% 
% Post Processing

    V_node_sq = zeros(N, 1);
    busesCounted = 0;
    
    %Sorting data based on actual bus numbers.
    for Area = 1:numAreas
        N_Area = N_Areas(Area);
        fprintf("Number of nodes in Area %d = %d.\n", Area, N_Area);
    
        if isRoot(Area) %Parent-less area will contribute all of its buses to the final tally
            v_area1 = v2(1:N_Area, Area);
            busesAdded = N_Area;
        else %Every child area will NOT double count the two buses it shares with its parent.
            v_area1 = v2(3:N_Area, Area);
            busesAdded = N_Area - 2;
        end
        
        
        V_node_sq(busesCounted+1: busesCounted + busesAdded) = v_area1;
        busesCounted = busesCounted + busesAdded;
    
    end
    
    V_node_sq_actual = zeros(N, 1);
    V_node_sq_actual( actualBusNums ) = V_node_sq;
    v1_actual = transpose(sqrt(V_node_sq_actual));
    
    microIterationLosses = microIterationLosses(1:itr, :);
    time_dist = time_dist(1:itr, :);
    
    macroIterationLossesTotal = kVA_B * sum(microIterationLosses, 2); %sum over each row
    %representing each iteration
    dist_timeTotal = sum(time_dist, 2);
    
    if displaySimulationResultPlots
        displayAndSaveSimulationResultPlots(itr, macroIterationLossesTotal, dist_timeTotal, R_max, v1_actual, saveSimulationResultPlots, Area, systemName, numAreas);
    end
    
    lineLoss_kW = macroIterationLossesTotal(itr);
    substationPower_kW = real(S_child(1))*kVA_B;
    timeToSolveOPFs_s = max(sum(time_dist));
    
    disp('------------------------------------------------------------')
    disp(['Line Loss: ', num2str(lineLoss_kW),' kW'])                       
    disp(['Substation Power: ', num2str(substationPower_kW),' kW'])
    disp(['Number of Iterations: ', num2str(itr)])
    disp(['Time to Solve: ', num2str(timeToSolveOPFs_s), 's'])
    disp('------------------------------------------------------------')
    
    programRunTime = toc(start);
    
    if saveSimulationResults
        saveResults(systemName, numAreas, lineLoss_kW, substationPower_kW, timeToSolveOPFs_s, itr, programRunTime);
    end
    
end