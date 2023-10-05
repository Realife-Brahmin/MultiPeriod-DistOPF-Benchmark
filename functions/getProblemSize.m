function [nLinEqnsT, nNonLinEqnsT, nVarsT] = getProblemSize(areaInfo, T)
    %GETPROBLEMSIZE Computes the sizes of linear equations, nonlinear equations, and variables.
    %
    % Syntax:
    %   [nLinEqnsT, nNonLinEqnsT, nVarsT] = getProblemSize(areaInfo, T)
    %
    % Description:
    %   This function calculates the total number of linear equations, nonlinear 
    %   equations, and optimization variables across a given time horizon T for an 
    %   optimization problem. It uses the provided areaInfo structure to obtain 
    %   details about the optimization model.
    %
    % Inputs:
    %   - areaInfo: Structure containing information about the optimization area.
    %               It should have fields:
    %                   * N_Area - Number of nodes in the area.
    %                   * m_Area - Number of branches in the area.
    %                   * nDER_Area - Number of DERs in the area.
    %                   * nBatt_Area - Number of batteries in the area.
    %   - T: Time horizon (integer) over which the problem is defined.
    %
    % Outputs:
    %   - nLinEqnsT: Total number of linear equations across the time horizon.
    %   - nNonLinEqnsT: Total number of nonlinear equations across the time horizon.
    %   - nVarsT: Total number of optimization variables across the time horizon.
    %

    % Extract information from areaInfo
    N_Area = areaInfo.N_Area;
    m_Area = areaInfo.m_Area;
    nDER_Area = areaInfo.nDER_Area;
    nBatt_Area = areaInfo.nBatt_Area;

    % Compute the number of variables for a single time step
    listNumVars1 = [m_Area*ones(1, 3), N_Area, nDER_Area, nBatt_Area*ones(1, 4)];
    nVars1 = sum(listNumVars1);
    nVarsT = nVars1 * T;

    % Compute the number of linear equations for a single time step
    listNumLinEqns1 = [m_Area*ones(1, 2), N_Area, nBatt_Area];
    nLinEqns1 = sum(listNumLinEqns1);
    nLinEqnsT = nLinEqns1 * T;

    % Compute the number of nonlinear equations for a single time step
    listNumNonlinEqns1 = m_Area;
    nNonLinEqns1 = sum(listNumNonlinEqns1);
    nNonLinEqnsT = nNonLinEqns1 * T;

end
