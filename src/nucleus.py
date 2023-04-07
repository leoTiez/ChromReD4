from typing import Tuple, Union, List, Iterable
import numpy as np
from scipy.ndimage.filters import gaussian_filter

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap


class Nucleus:
    def __init__(
            self,
            diff_rate: float = .5,
            diff_update_rate: float = .1,
            k_on: float = .4,
            k_off: float = .1,
            size_nucl: Iterable[int] = (100, 100),
            chromatin_len: int = 3600,
            p_init: float = .15
    ):
        self.diff_rate = diff_rate
        self.diff_update_rate = diff_update_rate
        self.k_on = k_on
        self.k_off = k_off
        self.size = np.array(list(size_nucl))
        self.chromatin = None
        self.chromatin_map = np.zeros(self.size, dtype='bool')
        self._init_chromatin(chromatin_len)
        self.associated = np.zeros(chromatin_len, dtype='bool')
        self.repaired = np.zeros(chromatin_len, dtype='bool')
        if 0 < p_init < 1:
            self.protein = np.random.choice([0, 1], size=size_nucl, p=[1. - p_init, p_init])
        if p_init < 0:
            self.protein = np.zeros(self.size)
            xmin, xmax = self.size[0] // 2 - int(.1 * self.size[0]), self.size[0] // 2 + int(.1 * self.size[0])
            ymin, ymax = self.size[1] // 2 - int(.1 * self.size[1]), self.size[1] // 2 + int(.1 * self.size[1])
            self.protein[xmin:xmax, ymin:ymax] = 1

    def _init_chromatin(self, chromatin_len: int = 1000, direction_change: float = .6):
        def rectify_boundaries(cp):
            return cp % self.size

        self.chromatin_map = np.zeros(self.size, dtype='bool')
        curr_pos = np.array([np.random.choice(self.size[0]), np.random.choice(self.size[1])])
        self.chromatin = np.zeros((chromatin_len, 2), dtype='uint16')
        self.chromatin[0, :] = curr_pos
        update = np.array([-1, 0])
        for i in range(1, chromatin_len):
            # self avoiding
            self.chromatin_map[self.chromatin[:i, 0], self.chromatin[:i, 1]] = True
            tmp_pos = rectify_boundaries(curr_pos + update)
            if self.chromatin_map[tmp_pos[0], tmp_pos[1]] or np.random.random() < direction_change:
                update_mask = np.zeros_like(self.chromatin_map, dtype='bool')
                xmin, xmax = (curr_pos[0] - 1) % self.size[0], (curr_pos[0] + 2) % self.size[0]
                ymin, ymax = (curr_pos[1] - 1) % self.size[1], (curr_pos[1] + 2) % self.size[1]
                if xmin < xmax:
                    if ymin < ymax:
                        update_mask[xmin:xmax, ymin:ymax] = True
                    else:
                        update_mask[xmin:xmax, ymin:] = True
                        update_mask[xmin:xmax, :ymax] = True
                else:
                    if ymin < ymax:
                        update_mask[xmin:, ymin:ymax] = True
                        update_mask[:xmax, ymin:ymax] = True
                    else:
                        update_mask[xmin:, ymin:] = True
                        update_mask[xmin:, :ymax] = True
                        update_mask[:xmax, ymin:] = True
                        update_mask[:xmax, :ymax] = True

                if np.any(np.logical_and(update_mask, ~self.chromatin_map)):
                    update_x, update_y = np.where(np.logical_and(update_mask, ~self.chromatin_map))

                    update_idx = np.random.choice(len(update_x))
                    curr_pos = np.array([update_x[update_idx], update_y[update_idx]])
                    curr_pos = rectify_boundaries(curr_pos)
                else:
                    new_x, new_y = np.where(self.chromatin_map == 0)
                    new_idc = np.random.choice(len(new_x))
                    curr_pos = np.array([new_x[new_idc], new_y[new_idc]])
            else:
                curr_pos = tmp_pos

            self.chromatin[i, :] = curr_pos

    def _free_closest(self, coord: np.ndarray, radius: float):
        x_absent, y_absent = np.where(self.protein == 0)
        absent_coord = np.asarray([x_absent, y_absent]).T
        dist = np.abs(absent_coord - coord)
        dist_torus = np.abs(self.size - dist)
        # Determine region w/ more movement
        idc = np.where(np.sum(np.minimum(dist, dist_torus), axis=1) < radius)[0]
        if len(idc) == 0:
            idc = np.array([np.argmin(np.sum(np.minimum(dist, dist_torus), axis=1))])
        return idc

    def _diffuse(self, sigma: List[float] = [.5, .5]):
        present = np.sum(self.protein == 1)
        max_dist = np.linalg.norm(self.size)
        decrease_idx = np.random.choice(present, size=int(self.diff_update_rate * float(present)))

        # Determine diffusability
        p_ratio = -gaussian_filter(self.protein.astype('float'), sigma, mode='constant')[self.protein == 1]
        p_ratio -= p_ratio.min()
        p_ratio = p_ratio[decrease_idx] / np.max(p_ratio[decrease_idx])
        x_present, y_present = np.where(self.protein == 1)
        for x, y, p in zip(x_present[decrease_idx], y_present[decrease_idx], p_ratio):
            if self.protein[x, y] == 0:
                continue
            if np.sum(self.protein == 1) != present:
                raise ValueError('Diffusion failed. Report a bug')
            x_absent, y_absent = np.where(self.protein == 0)
            increase_idc = self._free_closest(np.asarray([x, y]), max_dist * self.diff_rate * p)
            increase_idx = np.random.choice(increase_idc, size=1)
            self.protein[x, y] -= 1
            self.protein[x_absent[increase_idx], y_absent[increase_idx]] += 1

    def _associate(self):
        associate_mask = np.logical_and(
            (self.protein > 0)[self.chromatin[:, 0], self.chromatin[:, 1]],
            ~self.associated
        )
        association_mask = np.logical_and(associate_mask, np.random.random(len(associate_mask)) < self.k_on)
        self.associated[association_mask] = True
        self.repaired[association_mask] = True
        self.protein[self.chromatin[association_mask, 0], self.chromatin[association_mask, 1]] -= 1

    def _dissociate(self):
        dissociate_mask = self.associated.copy()
        dissociate_mask = np.logical_and(dissociate_mask, np.random.random(len(dissociate_mask)) < self.k_off)
        self.associated[dissociate_mask] = False
        self.protein[self.chromatin[dissociate_mask, 0], self.chromatin[dissociate_mask, 1]] += 1

    def reset(self, k_on: float, k_off: float):
        self.repaired[:] = False
        for i, coord in enumerate(self.chromatin[self.associated]):
            self.associated[i] = False
            self.protein[coord[0], coord[1]] = +1
        self.k_off = k_off
        self.k_on = k_on

    def update(self):
        self._diffuse()
        self._associate()
        self._dissociate()

    def plot(
            self,
            ax_plane: plt.Axes,
            ax_string: plt.Axes,
            reshape_chromatin: Tuple[int, int],
            t: int,
            plot_chromatin: bool = True
    ):
        ax_plane.clear()
        norm = Normalize(vmin=-.1, vmax=1.1)
        im = ax_plane.imshow(self.protein, cmap='coolwarm', norm=norm, alpha=1.)
        if plot_chromatin:
            colours = np.array([[0., 0., 0., 0.], [0., 0., 0., .5], [1., 0., 0., .5]])
            chromatin_map = self.chromatin_map.copy().astype('int8')
            chromatin_map[self.chromatin[self.associated, 0], self.chromatin[self.associated, 1]] = 2
            ax_plane.imshow(chromatin_map, cmap=ListedColormap(colours))
        ax_plane.set_title('Protein: %d' % t)
        ax_string.clear()
        im_chromatin = ax_string.imshow(self.repaired.reshape(reshape_chromatin), cmap='RdYlGn', norm=norm)
        ax_string.set_title('Visited DNA: %d' % t)

        return im, im_chromatin



