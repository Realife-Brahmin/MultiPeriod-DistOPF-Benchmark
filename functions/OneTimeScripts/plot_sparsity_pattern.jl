using Plots

function plot_sparsity_pattern(Aeq, beq)
    # Dark theme settings
    theme(:dark)
    
    # Get indices and values for non-zero elements in Aeq
    iAeq, jAeq, vAeq = [], [], []
    for i in 1:size(Aeq, 1)
        for j in 1:size(Aeq, 2)
            if Aeq[i, j] != 0
                push!(iAeq, size(Aeq, 1) - i + 1)  # Invert the y-values
                push!(jAeq, j)
                push!(vAeq, Aeq[i, j])
            end
        end
    end
    
    # Get indices for non-zero elements in beq
    nzBeq = findall(!iszero, beq)
    nzBeq = size(beq, 1) .- nzBeq .+ 1  # Invert the y-values for beq as well

    # Plotting
    p1 = scatter(jAeq[vAeq .< 0], iAeq[vAeq .< 0], markersize = 5, markerstrokewidth = 0, color = :red, label = "< 0", legend = :outertopright, xlims = (1, size(Aeq, 2)), title = "Sparsity pattern of Aeq", grid = true, xaxis = "Columns", yaxis = "Rows")
    scatter!(p1, jAeq[vAeq .== 0], iAeq[vAeq .== 0], markersize = 3, markerstrokewidth = 0, color = :white, label = "= 0")
    scatter!(p1, jAeq[vAeq .> 0], iAeq[vAeq .> 0], markersize = 5, markerstrokewidth = 0, color = :green, label = "> 0")

    p2 = scatter(repeat([1.5], length(nzBeq)), nzBeq, markersize = 5, markerstrokewidth = 0, color = :green, label = "beq", xlims = (1, 2), title = "Sparsity pattern of beq", grid = true, xaxis = " ", yaxis = "Rows", legend = :outertopright)

    # Combine plots side by side
    plot(p1, p2, layout = (1, 2))
end

A = [1 0 0; 0 1 0; 0 0 1]
b = [1, 0, 1]
plot_sparsity_pattern(A, b)
