function [c, ceq] = NonLinEqualities(x, Area, m_Area, fbus_Area, tbus_Area, indices_Pij, indices_Qij, indices_lij, indices_vAllj, itr, systemName, numAreas, varargin)
    
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
    
    strArea = convert2doubleDigits(Area);
    saveLocationFilename = strcat(saveLocation, systemName, "/numAreas_", num2str(numAreas), "/eqcons_area", strArea, fileExtension);
    fileOpenedFlag = false;

    if verbose && saveToFile && itr == 0 && Area == 2
        fileOpenedFlag = true;
        fid = fopen(saveLocationFilename, 'w');  % Open file for writing
    else
        verbose = false;
        fid = 1;
    end
    
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
        i_Idx = find(tbus_Area == j);
        i = fbus_Area(i_Idx);
    
        if ~isempty(i_Idx)
            myfprintf(verbose, fid, "It's parent is bus %d at index %d.\n", i, i_Idx);
            js_indices = find(fbus_Area == i);
            js = tbus_Area(js_indices);
            jes_Idx = js_indices(1);
            myfprintf(verbose, fid,  "Sibling(s) is(are):\n");
            myfprintf(verbose, fid, "%d ", js);
            myfprintf(verbose, fid, "\nlocated at:\n");
            myfprintf(verbose, fid, "%d ", js_indices);
            myfprintf(verbose, fid, "\nceq(%d) = l(%d) * v(%d) -  ( P(%d)^2 +  Q(%d)^2 )\n", i_Idx, i_Idx, js_indices(1), i_Idx, i_Idx)
            % myfprintf(verbose, fid, "ceq(%d) = x(%d) * x(%d) -  ( x(%d)^2 +  x(%d)^2 )\n", i_Idx, indices_l_Area(i_Idx), indices_v_Full_Area(jes_Idx), indices_P_Area(i_Idx), indices_Q_Area(i_Idx))
            
            indices_Pij_T = getIndicesT(indices_Pij, i_Idx);
            indices_Qij_T = getIndicesT(indices_Qij, i_Idx);
            indices_lij_T = getIndicesT(indices_lij, i_Idx);
            indices_vi_T = getIndicesT(indices_vAllj, jes_Idx);
            
            row = getIndices(indices_currentMag, i_Idx);
            
            ceq(row) = x(indices_lij_T) * x(indices_vi_T) -  ( x(indices_Pij_T)^2 +  x(indices_Qij_T)^2 );
        else
            myfprintf(verbose, fid, "It has NO parent bus.\n");
        end
    end

    if fileOpenedFlag
        fclose(fid);
    end

end