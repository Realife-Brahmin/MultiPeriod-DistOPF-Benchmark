function [c, ceq] = eqcons(x, Area, N_Area, fbus_Area, tbus_Area, indices_P_Area, indices_Q_Area, indices_l_Area, indices_v_Full_Area, itr, systemName, numAreas, varargin)
    
    verbose = false;
    saveToFile = false;
    saveLocation = "logfiles/";
    systemName = "ieee123";
    fileExtension = ".txt";

    for i = 1:2:numel(varargin)
        name = varargin{i};
        value = varargin{i+1};
        
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
    
    c = [];
    ceq = zeros(N_Area, 1);
    myfprintf(verbose, fid, "**********" + ...
        "Constructing ceq for Area %d.\n" + ...
        "***********\n", Area);
    
    for currentBusNum = 2 : N_Area
        myfprintf(verbose, fid, "*****\n" + ...
            "Checking for bus %d.\n" + ...
            "*****\n", currentBusNum);
        parentBusIdx = find(tbus_Area == currentBusNum);
        parentBusNum = fbus_Area(parentBusIdx);
    
        if ~isempty(parentBusIdx)
            myfprintf(verbose, fid, "It's parent is bus %d at index %d.\n", parentBusNum, parentBusIdx);
            siblingBusesIndices = find(fbus_Area == parentBusNum);
            siblingBuses = tbus_Area(siblingBusesIndices);
            eldestSiblingIdx = siblingBusesIndices(1);
            myfprintf(verbose, fid,  "Sibling(s) is(are):\n");
            myfprintf(verbose, fid, "%d ", siblingBuses);
            myfprintf(verbose, fid, "\nlocated at:\n");
            myfprintf(verbose, fid, "%d ", siblingBusesIndices);
            myfprintf(verbose, fid, "\nceq(%d) = l(%d) * v(%d) -  ( P(%d)^2 +  Q(%d)^2 )\n", parentBusIdx, parentBusIdx, siblingBusesIndices(1), parentBusIdx, parentBusIdx)
            myfprintf(verbose, fid, "ceq(%d) = x(%d) * x(%d) -  ( x(%d)^2 +  x(%d)^2 )\n", parentBusIdx, indices_l_Area(parentBusIdx), indices_v_Full_Area(eldestSiblingIdx), indices_P_Area(parentBusIdx), indices_Q_Area(parentBusIdx))
            % ceq(parentBusIdx) = x( indices_l_Area(parentBusIdx) ) * x( indices_v_Full_Area(eldestSiblingIdx) ) -  ( x( indices_P_Area(parentBusIdx) )^2 +  x( indices_Q_Area(parentBusIdx) )^2 );
            ceq(parentBusIdx) = x( indices_l_Area(parentBusIdx) ) * x( indices_v_Full_Area(parentBusNum) ) -  ( x( indices_P_Area(parentBusIdx) )^2 +  x( indices_Q_Area(parentBusIdx) )^2 );

        else
            myfprintf(verbose, fid, "It has NO parent bus.\n");
        end
    end

    if fileOpenedFlag
        fclose(fid);
    end

end