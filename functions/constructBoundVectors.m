% constructBoundVectors - Constructs two column vectors `lb` and `ub`
% containing lower and upper bounds for the specified power system variables
% based on the input arrays `arr`, `lbValues`, and `ubValues`.
%
% Syntax:
%   [lb, ub] = constructBoundVectors(arr, lbValues, ubValues)
%
% Input arguments:
% - arr: An array of t values representing the number of equations for
% different power system variables. For example, `arr = [1, m_Area-1, m_Area, m_Area, N_Area]`,
% where `m_Area` and `N_Area` are integers representing the number of equations for different variables.
%
% - lbValues: An array containing the minimum values for the power system variables.
% Each element of `lbValues` corresponds to a set of variables and their bounds.
% For example, `lbValues = {-1500, -1500, -1500, V_min^2}` represents the minimum values for `P`, `Q`, `l`, and `v`.
%
% - ubValues: An cell array containing the maximum values for the power system variables.
% Each element of `ubValues` corresponds to a set of variables and their bounds.
% For example, `ubValues = {1500, 1500, 1500, V_max^2}` represents the maximum values for `P`, `Q`, `l`, and `v`.
%
% Output arguments:
% - lb: A column vector representing the lower bounds for the power system variables.
%
% - ub: A column vector representing the upper bounds for the power system variables.
%
% Example:
%   N_Area = 19;
%   m_Area = 18;
%   arr = [1, m_Area-1, m_Area, m_Area, N_Area];
%   lbValues = {-1500, -1500, 0, V_min^2};
%   ubValues = {1500, 1500, 1500, V_max^2};
%   [lb, ub] = constructBoundVectors(arr, lbValues, ubValues);
%
% See also: (Add any relevant functions here if applicable)
function [lb, ub] = constructBoundVectors(arr, lbValues, ubValues)
    numLinOptEquations = sum(arr);
    if length(arr) ~= length(lbValues) || length(arr) ~= length(ubValues)
        error("Incompatible sizes of input.")
    end
    [lb, ub] = deal(zeros(numLinOptEquations, 1));

    indices = [0; transpose(cumsum(arr))];

    for i = 1:length(lbValues)
        lb(indices(i) + 1 : indices(i+1)) = lbValues(i);
        ub(indices(i) + 1 : indices(i+1)) = ubValues(i);
    end
end


