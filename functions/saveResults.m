function saveResults(systemName, numAreas, lineLoss_kW, substationPower_kW, timeToSolveOPFs_s, timePeriodNum, numMacroIterations, dopfRunTime)
    delta_t = 0.25;
    machineName = getenv("COMPUTERNAME");
    filenameSavedResults = strcat("processedData", filesep, systemName, filesep, "numAreas_", num2str(numAreas), filesep, "output_", machineName, ".txt");
    if timePeriodNum == 1
        fid = fopen(filenameSavedResults, 'w');
    else
        fid = fopen(filenameSavedResults, 'a');
    end
    
    fprintf(fid, strcat("Simulations run on machine: ", machineName));
    fprintf(fid, "Time Period = %d or Time = %.2f hours.\n", timePeriodNum, timePeriodNum*delta_t);
    fprintf(fid, "Line Loss = %.2fkW\n", lineLoss_kW);        
    fprintf(fid, "Substation Power = %.2fkW\n", substationPower_kW); 
    fprintf(fid, "Time to Solve = %.3fs\n", timeToSolveOPFs_s);    
    fprintf(fid, "Number of Macro-iterations = %d\n", numMacroIterations);
    fprintf(fid, "Run Time for this D-OPF: %.2f\n", dopfRunTime);
    fclose(fid);
end