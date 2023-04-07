import numpy as np
import matplotlib.pyplot as plt
import imageio
from pathlib import Path

from src.nucleus import Nucleus


def get_img(fig):
    graph = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    graph = graph.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    return graph


def main():
    n_iter_init = 300
    reshape = (60, 60)
    save_fig = True
    save_prefix = 'reaction_diffusion_koff_same'
    plot_chromatin = True
    nucl = Nucleus()
    plt.ion()
    fig_nucl, ax_nucl = plt.subplots(1, 1)
    fig_dna, ax_dna = plt.subplots(1, 1)
    nucl_fig_l, dna_fig_l = [], []
    for t in range(n_iter_init):
        if t == 150:
            nucl.reset(.4, .1)
        nucl.update()

        im_nucl, im_dna = nucl.plot(ax_nucl, ax_dna, reshape_chromatin=reshape, t=t, plot_chromatin=plot_chromatin)
        if t == 0:
            fig_nucl.colorbar(im_nucl)
            fig_dna.colorbar(im_dna)
        fig_nucl.tight_layout()
        fig_nucl.canvas.draw()
        fig_nucl.canvas.flush_events()
        fig_dna.tight_layout()
        fig_dna.canvas.draw()
        fig_dna.canvas.flush_events()
        nucl_fig_l.append(get_img(fig_nucl))
        dna_fig_l.append(get_img(fig_dna))

    plt.ioff()
    if save_fig:
        Path('figures/simulations').mkdir(exist_ok=True, parents=True)
        imageio.mimsave('figures/simulations/%s_nucleus.gif' % save_prefix, nucl_fig_l, fps=10)
        imageio.mimsave('figures/simulations/%s_dna.gif' % save_prefix, dna_fig_l, fps=10)
        plt.close('all')


if __name__ == '__main__':
    main()


