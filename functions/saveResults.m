function saveResults(systemName, numAreas, lineLoss_kW, substationPower_kW, timeToSolveOPFs_s, numIterations, programRunTime)
    filenameSavedResults = strcat("processedData/", systemName, "/numAreas_", num2str(numAreas), "/output.txt");
    fid = fopen(filenameSavedResults, 'w');
    fprintf(fid, "Line Loss = %.2fkW\n", lineLoss_kW);        
    fprintf(fid, "Substation Power = %.2fkW\n", substationPower_kW); 
    fprintf(fid, "Time to Solve = %.3fs\n", timeToSolveOPFs_s);    
    fprintf(fid, "Number of Iterations = %d\n", numIterations);
    fprintf(fid, "Program Run Time: %.2f", programRunTime);
    fclose(fid);
end