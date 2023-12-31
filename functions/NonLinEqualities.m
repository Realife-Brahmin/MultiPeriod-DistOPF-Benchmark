function [c, ceq] = NonLinEqualities(x, simInfo, areaInfo, T, varargin)
    
    noBatteries = simInfo.noBatteries;
    verbose = false;
    saveToFile = false;
    saveLocation = "logfiles/";
    systemName = "ieee123";
    fileExtension = ".txt";

    for kwarg_num = 1:2:numel(varargin)
        name = varargin{kwarg_num};
        value = varargin{kwarg_num+1};
        
        switch lower(name)
            case 'verbose'
                verbose = value;
            case 'savetofile'
                saveToFile = value;
            case 'savelocation'
                saveLocation = value;
            case 'systemname'
                systemName = value;
            case 'fileextension'
                fileExtension = value;
            otherwise
                error('Unknown option: %s', name);
        end
    end
    
    % strArea = convert2doubleDigits(areaInfo.Area);
    % saveLocationFilename = strcat(saveLocation, systemName, "/numAreas_", num2str(numAreas), "/eqcons_area", strArea, fileExtension);
    % fileOpenedFlag = false;

    if verbose && saveToFile && itr == 0 && Area == 2
        fileOpenedFlag = true;
        fid = fopen(saveLocationFilename, 'w');  % Open file for writing
    else
        verbose = false;
        fid = 1;
    end
    
    % Unpacking areaInfo
    N_Area = areaInfo.N_Area;
    m_Area = areaInfo.m_Area;
    nDER_Area = areaInfo.nDER_Area;
    if ~noBatteries
        nBatt_Area = areaInfo.nBatt_Area;
    else
        nBatt_Area = 0;
    end
    
    fb_Area = areaInfo.fb_Area;
    tb_Area = areaInfo.tb_Area;
    Area = areaInfo.Area;

    listNumVars1 = [m_Area*ones(1, 3), N_Area, nDER_Area, nBatt_Area*ones(1, 4)];

    varIndicesT = generateRangesFromValuesT(listNumVars1, T);
    

    indices_Pij = varIndicesT{1};
    indices_Qij = varIndicesT{2};
    indices_lij = varIndicesT{3};
    indices_vAllj = varIndicesT{4};

    listNumNonlinEqns1 = m_Area;
    nNonLinEqns1 = sum(listNumNonlinEqns1);
    nNonLinEqnsT = nNonLinEqns1 * T;
    eqnNonLinIndicesT = generateRangesFromValuesT(listNumNonlinEqns1, T);
    indices_currentMag = eqnNonLinIndicesT{1};
    
    ceq = zeros(nNonLinEqnsT, 1);
    c = [];
    % ceq = zeros(m_Area*T, 1);
    myfprintf(verbose, fid, "**********" + ...
        "Constructing ceq for Area %d.\n" + ...
        "***********\n", Area);
    
    for j = 2 : N_Area
        myfprintf(verbose, fid, "*****\n" + ...
            "Checking for bus %d.\n" + ...
            "*****\n", j);
        i_Idx = find(tb_Area == j);
        i = fb_Area(i_Idx);
    
        if ~isempty(i_Idx)
            % myfprintf(verbose, fid, "It's parent is bus %d at index %d.\n", i, i_Idx);
            js_indices = find(fb_Area == i);
            % js = tb_Area(js_indices);
            jes_Idx = js_indices(1);
            % myfprintf(verbose, fid,  "Sibling(s) is(are):\n");
            % myfprintf(verbose, fid, "%d ", js);
            % myfprintf(verbose, fid, "\nlocated at:\n");
            % myfprintf(verbose, fid, "%d ", js_indices);
            % myfprintf(verbose, fid, "\nceq(%d) = l(%d) * v(%d) -  ( P(%d)^2 +  Q(%d)^2 )\n", i_Idx, i_Idx, js_indices(1), i_Idx, i_Idx)
            % myfprintf(verbose, fid, "ceq(%d) = x(%d) * x(%d) -  ( x(%d)^2 +  x(%d)^2 )\n", i_Idx, indices_l_Area(i_Idx), indices_v_Full_Area(jes_Idx), indices_P_Area(i_Idx), indices_Q_Area(i_Idx))
            
            indices_Pij_T = getIndicesT(indices_Pij, i_Idx);
            indices_Qij_T = getIndicesT(indices_Qij, i_Idx);
            indices_lij_T = getIndicesT(indices_lij, i_Idx);
            % indices_vi_T = getIndicesT(indices_vAllj, jes_Idx);
            indices_vi_T = getIndicesT(indices_vAllj, i);

            
            row = getIndicesT(indices_currentMag, i_Idx);
            
            ceq(row) = x(indices_lij_T) .* x(indices_vi_T) -  ( x(indices_Pij_T).^2 +  x(indices_Qij_T).^2 );
        else
            myfprintf(verbose, fid, "It has NO parent bus.\n");
        end
    end

    % if fileOpenedFlag
    %     fclose(fid);
    % end

end