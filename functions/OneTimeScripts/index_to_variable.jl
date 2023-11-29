function index_to_variable(indices, areaInfo)
    m_Area = areaInfo[:m_Area]
    N_Area = areaInfo[:N_Area]
    nD_Area = areaInfo[:nD_Area]
    numVars = length(indices)
    notations = fill("", numVars)
    for varNum = 1:numVars
        index = indices[varNum]
        if index <= m_Area
            notations[varNum] = "Pij($(index))";
        elseif index <= 2 * m_Area
            notations[varNum] = "Qij($(index - m_Area))";
        elseif index <= 3 * m_Area
            notations[varNum] = "lij($(index - 2 * m_Area))";
        elseif index <= 3 * m_Area + N_Area
            notations[varNum] = "vj($(index - 3 * m_Area))";
        else
            notations[varNum] = "qD($(index - 3 * m_Area - N_Area))"
        end
    end

    return notations
end

# # Example usage:
# m_Area = 127
# N_Area = 128
# nD_Area = 85
# areaInfo = Dict(:m_Area => m_Area, :N_Area => N_Area, :nD_Area => nD_Area)
# index = [21, 275]  # Replace with the actual index
# variable_representation = index_to_variable(index, areaInfo)
# println(variable_representation)
