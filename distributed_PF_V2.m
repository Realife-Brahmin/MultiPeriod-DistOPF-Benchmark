
clear all
clc

%% Initialization:

loss_min = 1;   % CVR  = 0, Loss minimization = 1;

if (loss_min)
    CVR = [0 ,0];                       
else
    CVR = [0.6, 3.0];
end

% Connection Bus-
%col1: From Area, col2: To Area, col3: Connected dummy bus number%%

load CB.txt

Tot_Area = max(CB(:,2));

for k = 1:Tot_Area
    row(k) = size(find(k== CB(:,1)),1); %contains which area has how many children
end
%Check if the area is located at the end of the tree: 1- Yes, 0- No
end_area = ~row ;

%%
% Uncomment the line below for 13/123/smaller systems -
Se(1:max(CB(:,3))-1,1:Tot_Area) = 0;

Ss = zeros (Tot_Area,1);
S_A = zeros (max(row),Tot_Area);
Vs(1:Tot_Area) = 1.03^2;
Vs(1) = 1.03^2;
V(2,Tot_Area) = 0;
Dec_Var(1,Tot_Area) = 0;
%% Loop:

c=1;
iter = 0;
%%
while (c)
    %% Solve Areas parallelly:
    for Area = 1:Tot_Area
        [v2, Ss(Area), S_all, mic_iter(iter+1,Area), Dec_Var_1, time_dist(iter+1,Area), R] = ...
            NL_OPF_dist(Vs(Area), Se(CB(find(Area==CB(:,1)),3)-1,Area), Area, end_area(Area), loss_min,CVR);
        
        if (size(v2,1)< size (V,1))  % Ensuring the dimensions will match to store the voltages
            v2(end+1:size(V,1)) = 0;
            Dec_Var_1(end+1:size(V,1)) = 0;
            
        elseif (size(v2,1)> size (V,1))
            V(end+1:size(v2,1),:) = 0;
            Dec_Var(end+1:size(v2,1),:) = 0;
        end
        
        if (size(S_all,1)< size (S_A,1))  % Ensuring the dimensions will match to store the S flow
            S_all(end+1:size(S_A,1)) = 0;
            
        elseif (size(S_all,1)> size (S_A,1))
            S_A(end+1:size(S_all,1),:) = 0;
        end
        
        Dec_Var(:, Area) = (Dec_Var_1);
        V(:, Area) = (v2);   % Storing all the square of the node voltages
        S_A(:, Area) = (S_all);   % Storing all the S flow of the Area
        
        disp('*************')
        disp(['Iteration: ', num2str(iter), '    Area: ', num2str(Area)])
    end
    %% Stopping Criterion:
    eps = 0.001 ; % Tolerance
    
    % Residual Matrix:
    c = 0;
    R_max(iter+1) = 0;
    for k1 = 1:size(CB,1)
        R = [
            (S_A((CB(k1,3)-1),CB(k1,1))) - (Ss(CB(k1,2)));
            V(CB(k1,3),CB(k1,1)) - Vs(CB(k1,2))];        
        if max(abs(R))> R_max(iter+1)
            R_max(iter+1) = max(abs(R));
        end
    end
    if (R_max(iter+1)> eps)  % Checking if any boundary converged or not.
        c = 1;
        a(1) = 0.25;
        a(2) = 0.05;
        %% Communication-
        for k2 = 1:size(CB,1)
            V_iter(iter+1,CB(k2,3),CB(k2,1)) = V(CB(k2,3),CB(k2,1)) ;
            Se_iter(iter+1,CB(k2,2)) = Ss(CB(k2,2));
            
            if iter>2
                Vs(CB(k2,2)) = a(2)*V_iter(iter-1,CB(k2,3),CB(k2,1))+...
                    a(1)*V_iter(iter,CB(k2,3),CB(k2,1))+...
                    (1-sum(a))*V(CB(k2,3),CB(k2,1)) ;
                Se((CB(k2,3)-1),CB(k2,1)) = a(2)*Se_iter(iter-1,CB(k2,2))+...
                    a(1)*Se_iter(iter,CB(k2,2))+...
                    (1-sum(a))*Ss(CB(k2,2));
            else
                Vs(CB(k2,2)) = V(CB(k2,3),CB(k2,1)) ;
                Se((CB(k2,3)-1),CB(k2,1)) = Ss(CB(k2,2));
            end
        end
        iter = iter+1;
    end
    if (iter>150)
        c = 0;
    end
    
end
%% node voltages- 
disp('00000000000000')
V_node_sq = Vs(1);
Dec_Var_overall = 0;
for j = 1:Tot_Area
    %% min(CB((find(j==CB(:,1))),3))-1 = starting point of dummy buses
    v_area1 = V(2:min(CB((find(j==CB(:,1))),3))-1,j);
    DV_area1 = Dec_Var(2:min(CB((find(j==CB(:,1))),3))-1,j);
    
    if (end_area(j))
        end_node = min(find(V(:,j)==0));
        if (end_node)
            v_area1 = V(2:end_node-1,j);
            DV_area1 = Dec_Var(2:end_node-1,j);
            
            disp('---')
        else
            v_area1 = V(2:end,j);
            DV_area1 = Dec_Var(2:end,j);
            disp('%%%%%%')
        end
        %
    end
    node_size(j)= size(v_area1,1);
    V_node_sq = [V_node_sq;v_area1];
    Dec_Var_overall = [Dec_Var_overall;DV_area1 ];
    
    
end

%%
bus_ser = [1:36,119,126,127,117,54:68,122,118,37:53,120,125,128,69:116,121,123,124]';
% 
% bus_ser = [1:36,119,126,127,117,118,54:68,122,125,37:53,120,128,69:116,121,123,124];

for j = 1:length(bus_ser)
    V_node_sq_1(bus_ser(j),1) = V_node_sq(j);
    Dec_Var_overall_1(bus_ser(j),1) = Dec_Var_overall(j);
end
%%
V_node = (sqrt(V_node_sq_1))';
Dec_Var = Dec_Var_overall_1';

disp('------------------------------------------------------------')
if (loss_min)
    disp(['Line Loss: ', num2str(sum(mic_iter(end,:))*1000),' kW'])                       
% else 
    disp(['Substation Power: ', num2str(real(S_A(1))*1000),' kW']) 
end

disp(['Time to Solve: ', num2str(sum(max(time_dist,[],2))), 'sec'])
disp('------------------------------------------------------------')


