function getIdxFromRow(row, areaInfo)
    m_Area = areaInfo[:m_Area]
    N_Area = areaInfo[:N_Area]
    nD_Area = areaInfo[:nD_Area]
    if row <= m_Area
        return row
    elseif row <= 2 * m_Area
        return row - m_Area
    elseif row <= 2 * m_Area + N_Area
        return row - 2 * m_Area
    else 
        return 1
    end
end

# # Example usage
# m_Area = 127
# N_Area = 128
# nD_Area = 85
# areaInfo = Dict("m_Area" => m_Area, "N_Area" => N_Area, "nD_Area" => nD_Area)
# row = rand(1:2*m_Area+N_Area)  # Example input
# branchIdx = getBranchIdxFromRow(row, areaInfo)
