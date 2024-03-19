function [c, ceq] = computeCeq(x, indices_NonLin, indices_Pflow, indices_lij, indices_vAll_jes_T, indices_Pij_T, indices_Qij_T)
    
    c = [];
    % Pre-allocate memory for ceq
    ceq = zeros(length(indices_Pflow_ij_T)+1);
    
    % Calculate each ceq value
    ceq = x(indices_lij_T) .* x(indices_vAll_jes_T) - ...
          (x(indices_Pij_T).^2 + x(indices_Qij_T).^2);

end
