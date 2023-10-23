using Plots

function plotSparsity(Aeq::Matrix, beq::Vector)
    # Define custom colors
    wineRed = RGB(0.6, 0.2, 0.2)
    sexyGreen = RGB(0, 0.6, 0.3)
    darkGreen = RGB(0, 0.5, 0)
    dracula_bg = RGB(40/255, 42/255, 54/255)
    dracula_fg = RGB(248/255, 248/255, 242/255)

    # Extract the non-zero indices from Aeq and beq
    nzAeq = findall(!iszero, Aeq)
    nzBeq = findall(!iszero, beq)

    iAeq = [t[1] for t in nzAeq]
    jAeq = [t[2] for t in nzAeq]
    vAeq = Aeq[nzAeq]

    # Modify ticks for Aeq
    step_size_aeq = max(1, floor(Int, size(Aeq, 1) / 10))
    xticks_aeq = 0.5:step_size_aeq:size(Aeq, 2) + 0.5
    yticks_aeq = 0.5:step_size_aeq:size(Aeq, 1) + 0.5

    p1 = scatter(jAeq[vAeq .< 0], iAeq[vAeq .< 0], color=sexyGreen, legend=false, ratio=1, xlims=(0, size(Aeq,2)+1), ylims=(0, size(Aeq,1)+1), markersize=4, title="Sparsity pattern of Aeq", grid=true, background_color=dracula_bg, gridcolor=dracula_fg, xticks=xticks_aeq, yticks=yticks_aeq, xrotation=90)
    scatter!(p1, jAeq[vAeq .> 0], iAeq[vAeq .> 0], color=wineRed, markersize=4)

    # Modify ticks for beq with a larger step size
    step_size_beq = max(1, floor(Int, length(beq) / 10))
    xticks_beq = 0.5:1:1.5
    yticks_beq = 1:step_size_beq:length(beq)
    p2 = scatter(repeat([1], length(nzBeq)), nzBeq, color=darkGreen, xlims=(0, 2), ylims=(0, size(Aeq,1)+1), markersize=4, legend=false, ratio=1, title="Sparsity pattern of beq", grid=true, background_color=dracula_bg, gridcolor=dracula_fg, xticks=xticks_beq, yticks=yticks_beq)

    plot(p1, p2, layout=(1,2))
end


A = [1 0 1; 0 1 0; 1 0 1]
b = [1, 0, 1]
plotSparsity(A, b)