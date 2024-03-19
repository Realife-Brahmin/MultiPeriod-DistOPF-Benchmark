function f = initWithBatteries(x, indices_Pc, indices_Pd, nBatt_Area, varargin)
    
    f = 0;
    alpha = 0.5;
    etta_C = 0.8;
    etta_D = 0.8;

    for i = 1:nBatt_Area
        Pc_Idx = indices_Pc(i);
        Pd_Idx = indices_Pd(i);
        
        f = f + alpha* ( x(Pc_Idx)*(1-etta_C) + x(Pd_Idx)*(1/etta_D - 1) );
    end

end