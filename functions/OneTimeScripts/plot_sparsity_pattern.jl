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

    # Plotting for Aeq
    p = scatter(jAeq[vAeq .< 0], iAeq[vAeq .< 0], markersize = 5, markerstrokewidth = 0, color = :red, label = "Aeq < 0", legend = :outertopright, xlims = (0, size(Aeq, 2) + 3), ylims = (1, size(Aeq, 1)), title = "Sparsity pattern of Aeq and beq", grid = true, xaxis = "Columns and beq", yaxis = "Rows", aspect_ratio = :auto)
    scatter!(p, jAeq[vAeq .> 0], iAeq[vAeq .> 0], markersize = 5, markerstrokewidth = 0, color = :green, label = "Aeq > 0")

    # Plotting for beq in the same plot
    scatter!(p, repeat([size(Aeq, 2) + 2], length(nzBeq)), nzBeq, markersize = 5, markerstrokewidth = 0, color = :blue, label = "beq")

    return p
end

A = [1 0 0; 0 1 0; 0 0 1]
b = [1, 0, 1]
plot_sparsity_pattern(A, b)
