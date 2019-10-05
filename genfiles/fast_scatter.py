import matplotlib
import numpy             as np
import matplotlib.pyplot as plt

from   mpl_toolkits.axes_grid1 import   make_axes_locatable

def fast_scatter(ax, xs, ys, values, mmin, mmax, step, markersize=0.1):
  ##  Digitzie.                                                                                                                                                                                                                        
  points       = step * np.floor(np.clip(values, a_min=mmin, a_max=mmax) / step)
  levels       = np.unique(points)

  cmap         = plt.get_cmap("viridis", len(levels))
  norm         = matplotlib.colors.Normalize(vmin=levels[0], vmax=levels[-1])

  colors       = cmap([1. * x / len(levels) for x in range(len(levels))])

  for i, level in enumerate(levels):
    isin       = (points == level)

    print('Plotting level {} of {} - {} targets at {}.'.format(i, len(levels), np.count_nonzero(isin), level))

    ax.plot(xs[isin], ys[isin], markersize=markersize, c=colors[i], lw=0, marker='.')

  divider      = make_axes_locatable(ax)

  cax          = divider.append_axes('right', size='2%', pad=0.2)                                                                                                                                                 
  cb           = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)

  return  0
