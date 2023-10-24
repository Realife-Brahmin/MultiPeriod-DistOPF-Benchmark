using Plots
using GR

function plot_sparsity_pattern(
    Aeq, beq; 
    savePlots::Bool=false,
    saveLocation::String="processedData/",
    filenames=["filenameundefined.png"])

    theme(:dark)
    
    # Get indices and values for non-zero elements in Aeq
    nzAeq = findall(!iszero, Aeq)
    iAeq = size(Aeq, 1) .- getindex.(nzAeq, 1) .+ 1  # Invert the y-values
    jAeq = getindex.(nzAeq, 2)
    vAeq = Aeq[nzAeq]

    # Get indices and values for non-zero elements in beq
    nzBeq = findall(!iszero, beq)
    vBeq = beq[nzBeq]
    nzBeq .= size(beq, 1) .- nzBeq .+ 1  # Invert the y-values for beq as well

    # Create the scatter plots for Aeq
    p = scatter(jAeq[vAeq .< 0], iAeq[vAeq .< 0], markershape = :square, markersize = 5, color = :red, label = "Aeq < 0", legend = :outertopright, title = "Sparsity pattern of Aeq and beq", xaxis = "Columns and beq", yaxis = "Rows", aspect_ratio = :auto, grid = false)  # Set grid to false initially
    scatter!(jAeq[vAeq .> 0], iAeq[vAeq .> 0], markershape = :square, markersize = 5, color = :green, label = "Aeq > 0")

    # Plotting for beq
    scatter!(repeat([size(Aeq, 2) + 2], length(nzBeq[vBeq .< 0])), nzBeq[vBeq .< 0], markershape = :square, markersize = 5, color = :blue, label = "beq < 0")
    scatter!(repeat([size(Aeq, 2) + 2], length(nzBeq[vBeq .> 0])), nzBeq[vBeq .> 0], markershape = :square, markersize = 5, color = :cyan, label = "beq > 0")

    # Adjusting the ticks
    xt = 0:max(floor(0.1*size(Aeq, 2)), 1):size(Aeq, 2) + 2
    yt = 0:max(floor(0.1*size(Aeq, 2)), 1):size(Aeq, 1)
    xaxis!(p, ticks = (xt, string.(xt)))
    yaxis!(p, ticks = (yt, string.(yt)))

    # Add manually set gridlines
    plot!(p,  size=(1920, 1080), grid = (:both, :on, 0.5, :white, 0.5), xticks = xt, yticks = yt)  # Adjust the grid style as per your preference

    if savePlots
        gr()
        for filename ∈ filenames
            if filename == "filenameundefined.png"
                Plots.savefig(saveLocation*"Aeq_and_beq.png")
            else
                Plots.savefig(filename)
            end
        end
    end

    return p
end

# A = [1 0 0; 0 1 0; 0 0 1]
# b = [1, -1, 1]
# plot_sparsity_pattern(A, b)
