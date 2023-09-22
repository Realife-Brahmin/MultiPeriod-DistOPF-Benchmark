function [Aeq, beq] = LinEqualities()

    listNumVars1 = [m_Area*ones(1, 3), N_Area, nDER_Area, nBatt_Area*ones(1, 4)];
    nVars1 = sum(listNumVars1);
    nVarsT = nVars1 * T;
    listNumEqns1 = [m_Area*ones(1, 2), N_Area, nBatt_Area];
    nEqns1 = sum(listNumEqns1);
    nEqnsT = nEqns1 * T;
    
    Aeq = zeros(nEqnsT, nVarsT);
    sz = size(Aeq);
    beq = zeros(nEqnsT, 1);

    eqnIndicesT = generateRangesFromValuesT(listNumEqns1, T);
    varIndicesT = generateRangesFromValuesT(listNumVars1, T);
    

    indices_Pflow = eqnIndicesT{1};
    indices_Qflow = eqnIndicesT{2};
    indices_KVL = eqnIndicesT{3};
    indices_SOC = eqnIndicesT{4};

    indices_Pij = varIndicesT{1};
    indices_Qij = varIndicesT{2};
    indices_lij = varIndicesT{3};
    indices_vAllj = varIndicesT{4};
    indices_vj = excludeFirstElement(indices_vAllj);
    indices_qDj = varIndicesT{5};
    indices_Bj = varIndicesT{6};
    indices_Pdj = varIndicesT{7};
    indices_Pcj = varIndicesT{8};
    indices_qBj = varIndicesT{9};

     for j = 2 : N_Area
        myfprintf(true, 1, "*****\n" + ...
            "Checking for bus %d.\n" + ...
            "*****\n", j);       
 
        % The row index showing the 'parent' bus of our currentBus:
        
        i_Idx = find(tb_Area == j);
        i = fb_Area(i_Idx);
        myfprintf(logging_Aeq_beq, fid_Aeq_beq, "The parent of bus %d is bus %d at index %d.\n", j, i, i_Idx);
        
        k_indices = find(fb_Area == j);
        % disp(k_indices);
        ks = tb_Area(k_indices);

        js_indices = find(fb_Area==i);
        js = tb_Area(js_indices);

        jes_Idx = js_indices(1);
        jes = js(1);

        indices_Pflow_ij_T = getIndicesT(indices_Pflow, i_Idx);
        indices_Qflow_ij_T = getIndicesT(indices_Qflow, i_Idx);
        indices_KVL_ij_T = getIndicesT(indices_KVL, i_Idx);
        % indices_SOC_j_T = getIndicesT(indices_SOC, i_Idx)
        
        indices_Pij_T = getIndicesT(indices_Pij, i_Idx);
        indices_Qij_T = getIndicesT(indices_Qij, i_Idx);
        indices_lij_T = getIndicesT(indices_lij, i_Idx);
        indices_vAllj_T = getIndicesT(indices_vAllj, i_Idx);
        indices_vj_T = getIndicesT(indices_vj, i_Idx);
        
        % PFlow equations in Aeq, beq | BFM variables only
        row = indices_Pflow_ij_T;
        % indices_Pij_T
        Aeq(sub2ind(sz, row, indices_Pij_T)) = 1;
        for k_num = 1:length(k_indices)
            k_Idx = k_indices(k_num);
            indices_Pjk_T = getIndicesT(indices_Pij, k_Idx);
            Aeq(sub2ind(sz, row, indices_Pjk_T)) = -1;
        end
        Aeq(sub2ind(sz, row, indices_lij_T)) = -R_Area_Matrix(i, j);
        beq(row) = lambdaVals.*P_L_Area(j);
        
        % PIdx = i_Idx;
        % Aeq_Full( PIdx, indices_P(i_Idx) ) = 1;
        % Aeq_Full( PIdx, indices_l(i_Idx) ) = -R_Area_Matrix( i, j );
        % Aeq_Full( PIdx, indices_v(i_Idx) ) = -0.5 * CVR_P * P_L_Area( j );

        % QFlow equations in Aeq, beq | BFM variables only
        row = indices_Qflow_ij_T;
        % indices_Qij_T
        % indices_Qij

        Aeq(sub2ind(sz, row, indices_Qij_T)) = 1;
        for k_num = 1:length(k_indices)
            k_Idx = k_indices(k_num);
            indices_Qjk_T = getIndicesT(indices_Qij, k_Idx);
            Aeq(sub2ind(sz, row, indices_Qjk_T)) = -1;
        end
        Aeq(sub2ind(sz, row, indices_lij_T)) = -X_Area_Matrix(i, j);
        beq(row) = lambdaVals.*Q_L_Area(j) - Q_C_Area(j);

        %Q equations
        % QIdx = PIdx + m_Area;
        % Aeq_Full( QIdx, indices_Q(i_Idx) ) = 1;
        % Aeq_Full( QIdx, indices_l(i_Idx) ) = -X_Area_Matrix( i, j );
        % Aeq_Full( QIdx, indices_v(i_Idx) ) = -0.5 * CVR_Q * Q_L_Area( j );

        
       % List of Row Indices showing the set of 'children' buses 'under' our currentBus:
        % k_indices = find(fb_Area == j);
        % if ~isempty(k_indices)
        %     Aeq_Full(PIdx, indices_P(k_indices) ) = -1;   % for P
        %     Aeq_Full(QIdx, indices_Q(k_indices) ) = -1;   % for Q
        % end
        
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, P(%d)) = 1.\n", PIdx, i_Idx);
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, l(%d)) = -r(%d, %d).\n", PIdx, i_Idx, i, j);
        % for k_num = 1:length(k_indices)
        %     myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, P(%d)) = -1\n", PIdx, k_indices(k_num));
        % end
        % if CVR_P
        %     myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, v(%d)) = -0.5 * CVR_P * P_L(%d).\n", PIdx, i_Idx, j);
        % end
        

        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, Q(%d)) = 1.\n", QIdx, i_Idx);
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, l(%d)) = -x(%d, %d).\n", QIdx, i_Idx, i, j);
        % for k_num = 1:length(k_indices)
        %     myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, Q(%d)) = -1\n", QIdx, k_indices(k_num));
        % end
        % if CVR_Q
        %     myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, v(%d)) = -0.5 * CVR_Q * Q_L(%d).\n", QIdx, i_Idx, j);
        % end
        

        % KVL equations in Aeq, beq | Full equation (contains only BFM
        % variables)
        row = indices_KVL_ij_T;
        Aeq(sub2ind(sz, row, indices_vj_T)) = 1;
        % js_indices = find(fb_Area == i);
        % jes_Idx = js_indices(1);
        % jes = tb_Area(jes_Idx);
        indices_vi_T = getIndicesT(indices_vAllj, jes_Idx);
        Aeq(sub2ind(sz, row, indices_vi_T)) = -1;
        Aeq(sub2ind(sz, row, indices_lij_T)) = R_Area_Matrix(i, j)^2 + X_Area_Matrix(i, j)^2;
        Aeq(sub2ind(sz, row, indices_Pij_T)) = -2*R_Area_Matrix(i, j);
        Aeq(sub2ind(sz, row, indices_Qij_T)) = -2*X_Area_Matrix(i, j);
        
        % V equations
        % vIdx = QIdx + m_Area;
        % Aeq_Full( vIdx, indices_v(i_Idx) ) = 1;
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, v(%d)) = 1\n", vIdx, i_Idx);

        %Return the rows with the list of 'children' buses of 'under' the PARENT of our currentBus:
        %our currentBus will obviously also be included in the list.
        % js_indices = find(fb_Area == i);
        % js = tb_Area(js_indices);

        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "The siblings of bus %d\n", j);
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "include these buses: %d\n", js)
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "at indices %d.\n", js_indices);
        % jes_Idx = js_indices(1);
        % jes = js(1);
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq,  "which makes bus %d at index %d as the eldest sibling.\n", jes, jes_Idx);
        % Aeq_Full( vIdx, indices_vAll( jes_Idx ) ) = -1;
        % % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, v_Full(%d)) = -1\n", vIdx, jes_Idx);
        % Aeq_Full( vIdx, indices_P(i_Idx) ) = 2 * R_Area_Matrix( i, j );
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, P(%d)) = 2*r(%d, %d).\n", vIdx, i_Idx, i, j);
        % Aeq_Full( vIdx, indices_Q(i_Idx) ) = 2 * X_Area_Matrix( i, j );
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, Q(%d)) = 2*x(%d, %d).\n", vIdx, i_Idx, i, j);
        % Aeq_Full( vIdx, indices_l(i_Idx) ) = ...
            % -R_Area_Matrix( i, j )^2 + ...
            % -X_Area_Matrix( i, j )^2 ;
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, l(%d)) = -r(%d, %d)^2 -x(%d, %d)^2.\n", vIdx, i_Idx, i, j, i, j);
        

        % beq_Full(PIdx) = ...
            % ( 1- 0.5 * CVR_P ) * ...
            % ( P_L_Area( j ) - P_der_Area( j ) );
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "beq(%d) = (1 - 0.5*CVR_P)*(P_L(%d) - P_der(%d))\n", PIdx, j, j);
    
        % beq_Full(QIdx) =  ...
            % ( 1- 0.5*CVR_Q ) * ...
            % ( Q_L_Area( j ) - Q_C_Area( j ) );
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "beq(%d) = (1 - 0.5*CVR_Q)*(Q_L(%d) - Q_C(%d))\n", QIdx, j, j);

    end
    
    % 'Substation' Bus KVL | BFM variables required only
    indices_KVL_12_T = getIndicesT(indices_KVL, N_Area);
    indices_vAll1_T = getIndicesT(indices_vAllj, 1);
    row = indices_KVL_12_T;
    Aeq(sub2ind(sz, row, indices_vAll1_T)) = 1;
    beq(row) = v_parent_Area;

    % substation voltage equation
    % vSubIdx = 3*m_Area + 1;
    % Aeq_Full( vSubIdx, indices_vAll(1) ) = 1;
    % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, v_Full(1)) = 1\n", vSubIdx);

    % beq_Full(vSubIdx) = v_parent_Area;
    % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "beq(%d) = %.3f\n", vSubIdx, v_parent_Area);
    
    % DER equation addition
    % Table_DER = zeros(nDER_Area, 5);
    
    for der_num = 1:nDER_Area
        j = busesWithDERs_Area(der_num);
        % i_Idx = find(tb_Area == j);
        indices_qDj_T = getIndicesT(indices_qDj, der_num);
        
        row = indices_Pflow_ij_T;
        beq(row) = beq(row) - transpose(P_der_Area(j).*pvCoeffVals);

        row = indices_Qflow_ij_T;
        Aeq(sub2ind(sz, row, indices_qDj_T)) = 1;
        
        % QIdx = i_Idx + m_Area;
        % qD_Idx = indices_qD(der_num);
        % Aeq_Full(QIdx, qD_Idx) = 1;
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, qD(%d)) = 1\n", QIdx, der_num);
        
        %setting other parameters for DGs:
        % Table_DER(der_num, 2) = qD_Idx;
        
        % slope kq definiton:
        % Table_DER(der_num, 3) = 2*ub_qD_onlyDERbuses_Area(der_num)/(V_max-V_min); % Qmax at Vmin, and vice versa
        
        % Q_ref, V_ref definition:
        % Table_DER(der_num, 4) = Qref_DER;  %Qref
        % Table_DER(der_num, 5) = Vref_DER;  %Vref
    end
    
    for batt_num = 1:nBatt_Area

        j = busesWithBatts_Area(batt_num);
        % i_Idx = find(tb_Area == j);
        
        indices_SOC_j_T = getIndicesT(indices_SOC, batt_num);
        indices_SOC_j_T_2toT = indices_SOC_j_T(2:T);
        indices_SOC_j_T_1 = indices_SOC_j_T(1);
        
        indices_Bj_T = getIndicesT(indices_Bj, batt_num);
        indices_Pdj_T = getIndicesT(indices_Pdj, batt_num);
        indices_Pcj_T = getIndicesT(indices_Pcj, batt_num);
        indices_qBj_T = getIndicesT(indices_qBj, batt_num);
        
        % Pflow equations | Battery Variables
        row = indices_Pflow_ij_T;
        Aeq(sub2ind(sz, row, indices_Pdj_T)) = 1;
        Aeq(sub2ind(sz, row, indices_Pcj_T)) = -1;
        
        % Qflow equations | Battery Variables
        row = indices_Qflow_ij_T;
        Aeq(sub2ind(sz, row, indices_qBj_T)) = 1;
        
        % SOC equations | Battery Variables | First Time Interval
        row = indices_SOC_j_T_1;
        Aeq(row, indices_Bj_T(1)) = -1;
        Aeq(row, indices_Pdj_T(1)) = delta_t*etta_C;
        Aeq(row, indices_Pcj_T(1)) = -delta_t/etta_D;
        beq(row) = -B0Vals_pu_Area(batt_num);

        % SOC equations | Battery Variables | First Time Interval
        row = indices_SOC_j_T_2toT;
        Aeq(sub2ind(sz, row, indices_Bj_T(2:T))) = -1;
        Aeq(sub2ind(sz, row, indices_Bj_T(1:T-1))) = 1;
        Aeq(sub2ind(sz, row, indices_Pcj_T(2:T))) = delta_t*etta_C;
        Aeq(sub2ind(sz, row, indices_Pdj_T(2:T))) = -delta_t/etta_D;

        % PEqnIdx = i_Idx;
        % QEqnIdx = i_Idx + m_Area;
        % BEqnIdx = numLinOptEquationsBFM + batt_num;
        % 
        % B_Idx = indices_B(batt_num);
        % Pc_Idx = indices_Pc(batt_num);
        % Pd_Idx = indices_Pd(batt_num);
        % qB_Idx = indices_qB(batt_num);

        % Aeq_Full(PEqnIdx, Pc_Idx) = -1;
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, Pc(%d)) = -1\n", PEqnIdx, batt_num);

        % Aeq_Full(PEqnIdx, Pd_Idx) = 1;
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, Pd(%d)) = 1\n", PEqnIdx, batt_num);

        % Aeq_Full(QEqnIdx, qB_Idx) = 1;
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, qB(%d)) = 1\n", QEqnIdx, batt_num);
        
        % Aeq_Full(BEqnIdx, B_Idx) = 1;
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, B(%d)) = 1\n", BEqnIdx, batt_num);
        % Aeq_Full(BEqnIdx, Pc_Idx) = -delta_t*etta_C;
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, Pc(%d)) = -delta_t*etta_C\n", BEqnIdx, batt_num);
        % Aeq_Full(BEqnIdx, Pd_Idx) = delta_t/etta_D;
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "Aeq(%d, Pd(%d)) = delta_t*etta_D\n", BEqnIdx, batt_num);

        % beq_Full(BEqnIdx) = B0Vals_pu_Area(batt_num);
        % myfprintf(logging_Aeq_beq, fid_Aeq_beq, "beq(%d) = B0(%d) = %f\n", BEqnIdx, batt_num, B0Vals_pu_Area(batt_num));
    end

    if fileOpenedFlag_Aeq_beq
        fclose(fid_Aeq_beq);
    end
    
    plotSparsity(Aeq, beq);
end