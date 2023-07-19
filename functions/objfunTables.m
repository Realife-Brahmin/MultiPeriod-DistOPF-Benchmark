% objfunTables: Computes the value for a loss minimization objective function for a particular area.
%
% Syntax:
%   f = objfunTables(x, N_Area, fbus_Area, tbus_Area, indices_l_Area, R_Area_Matrix)
%   Using @(x)objfunTables, the function can be used to output the
%   objective function in terms of branch flow variables x: f(x), so that x can
%   later be solved for in order to minimize f(x).
%
% Input Arguments:
%   - x: Vector of branch flow values.
%   - N_Area: Total number of buses in the area.
%   - fbus_Area: Vector of "from" buses for the branches in the area.
%   - tbus_Area: Vector of "to" buses for the branches in the area.
%   - indices_l_Area: Indices of branches in the area.
%   - R_Area_Matrix: Resistance matrix for the branches in the area.
%
% Output Argument:
%   - f: Value of the loss minimization objective function for the area.
%
% Description:
%   This function computes the loss minimization objective function for a particular area in the branch flow model. The objective function is calculated by summing the products of branch flow values and corresponding resistances for each branch in the area.
%   The variables used in this function have the following meanings:
%   - x: Represents the branch flow values.
%   - N_Area: Represents the total number of buses in the area.
%   - fbus_Area: Represents the "from" buses for the branches in the area.
%   - tbus_Area: Represents the "to" buses for the branches in the area.
%   - indices_l_Area: Represents the indices of branches in the area.
%   - R_Area_Matrix: Represents the resistance matrix for the branches in the area.
%
% Example:
%   N_Area = 5;
%   fbus_Area = [1, 2, 2, 3, 4];
%   tbus_Area = [2, 3, 4, 5, 5];
%   indices_l_Area = [1, 2, 3, 4, 5];
%   R_Area_Matrix = [0, 0, 0, 0, 0; 
%                    0, 0.2, 0.3, 0, 0;
%                    0, 0, 0.1, 0.4, 0;
%                    0, 0, 0, 0.5, 0.6;
%                    0, 0, 0, 0, 0.7];
%   x = [0.1, 0.2, 0.3, 0.4, 0.5];
%   f = objfunTables(x, N_Area, fbus_Area, tbus_Area, indices_l_Area, R_Area_Matrix);
%
%   The output 'f' will be 0.13.
%
% See also: (Any related functions or files that are useful to mention)
%
% References: (If applicable)

function f = objfunTables(x, N_Area, fbus_Area, tbus_Area, indices_l_Area, R_Area_Matrix)

f = 0;

for currentBusNum = 2 : N_Area
    parentBusIdx = find(tbus_Area == currentBusNum);
    parentBusNum = fbus_Area(parentBusIdx);
    
    if ~isempty(parentBusIdx)
       f = f + x( indices_l_Area(parentBusIdx) ) * R_Area_Matrix( parentBusNum, currentBusNum );
    end
    
end