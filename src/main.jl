using ChromReD4
import Plots; const plt = Plots;

function main()
    len_chromatin = 10000
    size_nucleus = [50, 50]
    upd_steps = 150

    diffu_const = .245
    asso_no = .5
    disso_no = .3
    asso_rad = .3
    disso_rad = .05

    save_gif = false
    show_plots = true
    frames = []
    nuc = init_nucleus(len_chromatin, size_nucleus)
    for _ in 1:upd_steps
        update!(nuc, asso_prob=asso_no, disso_prob=disso_no, diffu=diffu_const)
        graph = plot(nuc, "Before Radiation")
        push!(frames, graph)
        if show_plots
            plt.display(graph)
            plt.gui()
        end
        sleep(0.1)
    end

    for _ in 1:upd_steps
        update!(nuc, asso_prob=asso_rad, disso_prob=disso_rad, diffu=diffu_const)
        graph = plot(nuc, "After Radiation")
        push!(frames, graph)
        if show_plots
            plt.display(graph)
            plt.gui()
        end
        sleep(0.1)
    end

    if save_gif
        animation_gif = plt.Animation()
        for i in size(frames)
            plt.frame(animation_gif, frames[i])
        end
        plt.gif(animation_gif, "ChromReD4_example.gif", fps = 30)
    end
end

main()