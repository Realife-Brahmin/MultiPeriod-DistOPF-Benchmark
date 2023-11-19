function resMax = computeBiggestResidual(y_allR_vs_macroItr, sysInfo, simInfo)
    numRelationships = sysInfo.numRelationships;
    macroItr = simInfo.macroItr + 1; % starts at 1

    y_allR_mI = y_allR_vs_macroItr(macroItr);
    y_allR_mIm1 = y_allR_vs_macroItr(macroItr-1);

    residualMax_allR_mI = 0;

    for Rnum = 1:numRelationships
        y_R_mI = y_allR_mI{Rnum};
        y_R_mIm1 = y_allR_mIm1{Rnum};
        residuals_R_mI = y_R_mI - y_R_mIm1;
        residualMax_R_mI = max(max(residuals_R_mI));
        residualMax_allR_mI = max(residualMax_allR_mI, residualMax_R_mI);
    end
    
    resMax = residualMax_allR_mI;
end