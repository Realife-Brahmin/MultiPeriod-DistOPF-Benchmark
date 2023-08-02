function f = objfun(x, N_Area, nDER_Area, nBatt_Area, fb_Area, tb_Area, R_Area_Matrix, X_Area_Matrix, varargin)
    
 % Default values for optional arguments
    verbose = false;
    etta_C = 0.80;
    etta_D = 0.80;
    % alpha = 0.5;
    alpha = 3e-5;
    mainObjFun = "func_PLoss";
    secondObjFun = "func_SCD";
    m_Area = length(fb_Area);
    indices_P = 1:m_Area;
    indices_Q = indices_P + m_Area;
    indices_l = 2*m_Area + 1:3*m_Area;
    indices_vFull_NoLoss = indices_Q + N_Area;
    % indices_v_NoLoss = indices_vFull_NoLoss(2:end);
    indices_Pc_NoLoss = 2*m_Area + N_Area + nDER_Area + nBatt_Area + 1:2*m_Area + N_Area + nDER_Area + 2*nBatt_Area;
    indices_Pd_NoLoss = indices_Pc_NoLoss + nBatt_Area;

    indices_Pc_Full = indices_Pc_NoLoss + m_Area;
    indices_Pd_Full = indices_Pd_NoLoss + m_Area;
    % Process optional arguments
    numArgs = numel(varargin);

    if mod(numArgs, 2) ~= 0
        error('Optional arguments must be specified as name-value pairs.');
    end
    
    validArgs = ["verbose", "etta_C", "etta_D", "alpha", "mainObjFun", ...
        "secondObjFun", "indices_P", "indices_Q", "indices_vFull_NoLoss", "indices_l", "indices_Pc", "indices_Pd",...
        "voltageVals"];
    
    for i = 1:2:numArgs
        argName = varargin{i};
        argValue = varargin{i+1};
        
        if ~ischar(argName) || ~any(argName == validArgs)
            error('Invalid optional argument name.');
        end
        
        switch argName
            case "verbose"
                verbose = argValue;
            case "etta_C"
                etta_C = argValue;
            case "etta_D"
                etta_D = argValue;
            case "mainObjFun"
                mainObjFun = argValue;
            case "secondObjFun"
                secondObjFun = argValue;
            case "alpha"
                alpha = argValue;
            case "indices_P"
                indices_P = argValue;
            case "indcies_Q"
                indices_Q = argValue;
            case "indices_l"
                indices_l = argValue;
            case "indices_vFull_NoLoss"
                indices_vFull_NoLoss = argValue;
            case "indices_Pc"
                indices_Pc = argValue;
            case "indices_Pd"
                indices_Pd = argValue;
            case "voltageVals"
                v0Vals = argValue;
        end
    end
    

    f = 0;
    
    if strcmp(mainObjFun, "func_PLoss")
        for currentBusNum = 2 : N_Area
            parentBusIdx = find(tb_Area == currentBusNum);
            parentBusNum = fb_Area(parentBusIdx);
            
            if ~isempty(parentBusIdx)
               f = f + x( indices_l(parentBusIdx) ) * R_Area_Matrix( parentBusNum, currentBusNum );
            end
            
        end
    elseif strcmp(mainObjFun, "func_PLoss-fake")
          if ~exist("v0Vals", "var")
              error("Fake Loss Minimization to be performed, but v0Vals NOT inserted.")
          end
          P0 = x(indices_P);
          Q0 = x(indices_Q);
          v0 = v0Vals;    
          l0_fake = zeros(m_Area, 1);
          for currentBusNum = 2 : N_Area
              parentBusIdx = find(tb_Area == currentBusNum);
              parentBusNum = fb_Area(parentBusIdx);
              siblingBusesIndices = find(parentBusNum == fb_Area);
              l0_fake( parentBusIdx ) = ( P0(parentBusIdx)^2 + Q0(parentBusIdx)^2 ) / v0(siblingBusesIndices(1));

              f = f + l0_fake(parentBusIdx) * R_Area_Matrix( parentBusNum, currentBusNum );
          end
    elseif strcmp(mainObjFun, "func_QLoss")
        l = x(indices_l);
        for currentBusNum = 2:N_Area
            parentBusIdx = find(tb_Area == currentBusNum);
            parentBusNum = fb_Area(parentBusIdx);
            f = f + l(parentBusIdx) * X_Area_Matrix(parentBusNum, currentBusNum);
        end
    elseif strcmp(mainObjFun, "none")
        myfprintf(verbose, "Okay, no primary objective function component.\n");
    else
        error("Unknown Primary Objective Function");
    end
    
    if strcmp(secondObjFun, "func_SCD")
        myfprintf(verbose, "Complementarity Needs to Be Respected.\n");

        if ~exist('indices_Pc', 'var') && strcmp(mainObjFun, "func_PLoss")
            indices_Pc = indices_Pc_Full;
            indices_Pd = indices_Pd_Full;
        elseif ~exist('indices_Pc', 'var') && ~strcmp(mainObjFun, "func_PLoss")
            indices_Pc = indices_Pc_NoLoss;
            indices_Pd = indices_Pd_NoLoss;
        elseif strcmp(mainObjFun, "none")
            error("Illegal! Cannot solve for SCD with no minimization function. Abort!");
        else
            myfprintf(verbose, "Hope that Pc and Pd indices were input as optional arguments, cause this shit is being solved for.\n");
        end

        for i = 1:nBatt_Area
            Pc_Idx = indices_Pc(i);
            Pd_Idx = indices_Pd(i);
            
            f = f + alpha* ( x(Pc_Idx)*(1-etta_C) + x(Pd_Idx)*(1/etta_D - 1) );
        end
    elseif strcmp(secondObjFun, "none")
        myfprintf(verbose, "No additional objective function.\n");
    else
        error("Unknown Secondary Objective Function");
    end

end