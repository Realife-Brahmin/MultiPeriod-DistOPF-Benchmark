function [nVars_allT, nLinConstraints_allT,  nNonLinEqns_allT] = getProblemSizeQuick(N, m, nDER, nBatt, T, varargin)
    % getProblemSizeQuick Computes the problem size for optimization
    % Inputs:
    %   N, m, nDER, nBatt - Problem parameters
    %   T - Time horizon
    %   Optional:
    %     'batteryTerminalChargeConstraint' - Constraint type for battery terminal ('soft' or 'hard'), default 'hard'
    
    p = inputParser;
    addParameter(p, 'batteryTerminalChargeConstraint', 'hard', @(x) any(validatestring(x, {'soft', 'hard'})));
    parse(p, varargin{:});
    
    batteryTerminalChargeConstraint = p.Results.batteryTerminalChargeConstraint;
    
    % Calculate number of variables at each time step and for all T
    listNumVars_t_1toT = [m*ones(1, 3), N, nDER, nBatt*ones(1, 4)];
    nVars_t_1toT = sum(listNumVars_t_1toT);
    nVars_allT = nVars_t_1toT * T;
    
    % Calculate number of linear equations for time steps 1 to T-1
    listNumLinEqns_t_1toTm1 = [m*ones(1, 2), N, nBatt];
    nLinEqns_t_1toTm1 = sum(listNumLinEqns_t_1toTm1);
    

    % Handle the terminal time step based on battery terminal charge constraint
    if strcmp(batteryTerminalChargeConstraint, "soft")
        listNumLinEqns_t_T = [m*ones(1, 2), N, nBatt];
    else  % Default to "hard" if not "soft"
        listNumLinEqns_t_T = [m*ones(1, 2), N, 2*nBatt];
    end
    nLinEqns_t_T = sum(listNumLinEqns_t_T);
    
    % Total number of linear equations
    nLinEqns_allT = nLinEqns_t_1toTm1 * (T-1) + nLinEqns_t_T;
    
    listNumLinIneqns_lb_t_1toT = [1, 0, m, N, nDER, nBatt*ones(1, 4)]; % P, Q, l, v, qD, B, Pc, Pd, qB
    listNumLinIneqns_ub_t_1toT = [0, 0, m, N, nDER, nBatt*ones(1, 4)]; % P, Q, l, v, qD, B, Pc, Pd, qB
    nLinIneqns_t_1toT = sum(listNumLinIneqns_lb_t_1toT) + sum(listNumLinIneqns_ub_t_1toT);
    nLinIneqns_allT = T * nLinIneqns_t_1toT;

    nLinConstraints_allT = nLinEqns_allT + nLinIneqns_allT;
    % Assuming non-linear equations count
    nNonLinEqns_allT = m * T;  % only P^2 + Q^2 - l*v == 0
    
end