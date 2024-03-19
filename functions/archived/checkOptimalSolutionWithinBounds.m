function flag = checkOptimalSolutionWithinBounds(optimalSolution, lb, ub)
    tolerance = 1;
    % Check if the optimal solution is within the bounds with a given tolerance
    
    % Ensure optimalSolution, lb, and ub are column vectors
    optimalSolution = optimalSolution(:);
    lb = lb(:);
    ub = ub(:);
    
    % Check if optimalSolution is within tolerance percent of the bounds
    range = ub - lb;
    tolerance_value = tolerance * range / 100;
    
    flag = all(optimalSolution >= lb + tolerance_value) && all(optimalSolution <= ub - tolerance_value);
end
