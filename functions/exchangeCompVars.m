function [P_L_Parent, Q_L_Parent] = exchangeCompVars(P_L_Parent, Q_L_Parent, S_connection_Child)
    % EXCHANGECOMPVARS exchanges results from OPF iterations of the child area to the parent area.
    %
    % Input:
    %   - P_L_Parent: Pre-existing active power values in the parent area.
    %   - Q_L_Parent: Pre-existing reactive power values in the parent area.
    %   - S_connection_Child: Complex power from the child area for interconnection nodes.
    %
    % Output:
    %   - P_L_Parent: Updated active power values in the parent area after receiving data from the child.
    %   - Q_L_Parent: Updated reactive power values in the parent area after receiving data from the child.
    
    numChildAreas = size(S_connection_Child, 1); % Could even be zero for a child-less area
    
    for childArea_num = 1:numChildAreas
        % The load of the parent area at the node of interconnection is
        % basically the interconnection area power
        P_L_Parent(end-childArea_num+1) = real(S_connection_Child(end-childArea_num+1));   % in PU
        Q_L_Parent(end-childArea_num+1) = imag(S_connection_Child(end-childArea_num+1));
    end
end
