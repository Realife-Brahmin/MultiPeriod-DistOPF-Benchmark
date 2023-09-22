function areaInfo = exchangeCompVars(areaInfo, S_connection_Area)
    %EXCHANGECOMPVARS Updates the P_L_Area and Q_L_Area in areaInfo based on S_connection_Area.
    %
    % Description:
    %   Given the S_connection_Area which contains the power information of interconnection nodes,
    %   this function updates the respective loads P_L_Area and Q_L_Area in the areaInfo structure.
    %
    % Input:
    %   - areaInfo: Structure containing various power system parameters of the area.
    %   - S_connection_Area: Complex power values at the interconnection nodes.
    %
    % Output:
    %   - areaInfo: Updated areaInfo structure.
    %

    numChildAreas = size(S_connection_Area, 1); % Could even be zero for a child-less area

    for childArea_num = 1:numChildAreas
        % The load of the parent area at the node of interconnection is
        % basically the interconnection area power
        areaInfo.P_L_Area(end-childArea_num+1) = real(S_connection_Area(end-childArea_num+1));
        areaInfo.Q_L_Area(end-childArea_num+1) = imag(S_connection_Area(end-childArea_num+1));
    end
end
