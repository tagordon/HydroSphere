import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches
import cmasher as cmr
import numpy as np

default_cmap = cmr.get_sub_cmap('cmr.arctic_r', 0.0, 0.8)

def plotphase(fig_or_ax, phase_array, x, y, xlabel=None, ylabel=None, name=None, cmap=default_cmap):
    
    is_figure = isinstance(fig_or_ax, mpl.figure.Figure)
    if is_figure:
        ax = plt.gca()
        fig = fig_or_ax
    else:
        ax = fig_or_ax
    
    phase_array[phase_array > 4] = phase_array[phase_array > 4] - 1
    ax.pcolormesh(x, y, phase_array.T, cmap=cmap, vmin=0, vmax=6,rasterized=True)
    
    if is_figure:
        
        ax.set_ylabel(xlabel, fontsize=15)
        ax.set_xlabel(ylabel, fontsize=15)
        ax.set_title('Phase of water in contact with rocky core', fontsize=20)
        
        #colorbar stuff 
        bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        cax = fig.add_axes([0.15, 0.03, 0.71, 0.03])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')
        cb.set_ticks([0, 1, 2, 3, 4, 5, 6])
        cb.set_ticklabels(['liquid water', 'ice-Ih', 'ice-II', 'ice-III', 'ice-V', 'ice-VI', 'ice-VII'], fontsize=15)
        cb.ax.tick_params(size=0)
        
    if name is not None:
        
        ax.annotate(
            name, 
            (0.6, 0.85), 
            xycoords='axes fraction', 
            fontsize=20, 
            bbox=dict(facecolor=cmap(1.0), alpha=1.0, edgecolor=None, pad=10.0)
        ) 
        ax.annotate(
            name, 
            (0.6, 0.85), 
            xycoords='axes fraction', 
            fontsize=20, 
            bbox=dict(facecolor='white', alpha=0.5, edgecolor=None, pad=10.0)
        )
        
    return ax

def plotphase_multiple(fig, axs, phase_arrays, x, y, xlabel=None, ylabel=None, names=None, cmap=default_cmap):
    
    axs_flat = axs.flatten()
    
    for ax, phases, n in zip(axs_flat, phase_arrays, names):
        plotphase(ax, phases, x, y, xlabel=xlabel, ylabel=ylabel, name=n, cmap=cmap)
        
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cax = fig.add_axes([0.15, 0.03, 0.71, 0.03])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')
    cb.set_ticks([0, 1, 2, 3, 4, 5, 6])
    cb.set_ticklabels(['liquid water', 'ice-Ih', 'ice-II', 'ice-III', 'ice-V', 'ice-VI', 'ice-VII'], fontsize=15)
    cb.ax.tick_params(size=0)
    
    [
        ax.annotate(
            name, 
            (0.6, 0.85), 
            xycoords='axes fraction', 
            fontsize=20, 
            bbox=dict(facecolor=cmap(1.0), alpha=1.0, edgecolor=None, pad=10.0)
        ) 
        for name, ax in zip(names, axs_flat)
    ]


    [
        ax.annotate(
            name, 
            (0.6, 0.85), 
            xycoords='axes fraction', 
            fontsize=20, 
            bbox=dict(facecolor='white', alpha=0.5, edgecolor=None, pad=10.0)
        ) 
        for name, ax in zip(names, axs_flat)
    ]
    
    [ax.grid() for ax in axs_flat]
    [ax.set_xlim(0.001 - np.mean(np.diff(x))/2, 0.03) for ax in axs_flat]
    [ax.set_ylim(0.5, 1.5) for ax in axs_flat]
    
    [ax.set_ylabel(r'density of rocky core ($\rho_\oplus$)', fontsize=15) for ax in axs[:,0]]
    [ax.set_xlabel('water mass fraction', fontsize=15) for ax in axs[-1,:]]

    fig.suptitle('Phase of water in contact with rocky core', fontsize=20)
    fig.subplots_adjust(top=0.92)

    cax = fig.add_axes([0.15, 0.03, 0.71, 0.03])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')
    cb.set_ticks([0, 1, 2, 3, 4, 5, 6])
    cb.set_ticklabels(['liquid water', 'ice-Ih', 'ice-II', 'ice-III', 'ice-V', 'ice-VI', 'ice-VII'], fontsize=15)
    cb.ax.tick_params(size=0)
    
def plot_wedge(ax, profile, extent=30, cmap=default_cmap):
    
    z = np.linspace(0, 1, len(profile))
    edge_ids = np.where(np.abs(np.diff(profile)) > 0)
    edges = z[edge_ids]
    phases = np.concatenate([profile[edge_ids], [profile[-1]]])[::-1]
    
    wedges = [
        patches.Wedge(
            (0.5, 0), 
            0.5 + d / 2, 
            90 - extent, 
            90 + extent, 
            facecolor=cmap(p/7), 
            edgecolor='k', 
            linewidth=1
        ) for p, d in zip(phases, edges)
    ]
    
    wedges.append(
        patches.Wedge(
            (0.5, 0), 
            1, 
            90 - extent, 
            90 + extent, 
            facecolor=cmap(profile[0]/7), 
            edgecolor='k', 
            linewidth=1)
    )
    
    wedges = wedges[::-1]
    
    ices = [
        'liquid water', 
        'ice-Ih', 
        'ice-II', 
        'ice-III', 
        'ice-IV', 
        'ice-V', 
        'ice-VI', 
        'ice-VII'
    ]
    ylocs = np.concatenate([[0], edges]) + np.diff(np.concatenate([[0], edges, [1]]))/3
    ylocs = 0.5 + ylocs / 2
    xloc = lambda y: 0.55 + np.tan(extent * np.pi / 180) * y 
    ylocs = np.cos(extent * np.pi / 180) * ylocs
    
    [ax.annotate(ices[i], xy=(xloc(y), y), xycoords='data', fontsize=20) for i, y in zip(np.int64(phases), ylocs)]
    ax.annotate('rocky core', xy=(0.43, 0.4), xycoords='data', fontsize=20)

    ax.set_ylim(0.3, 1.1)
    ax.set_xlim(0.5 - np.tan(extent * np.pi / 180), 0.5 + np.tan(extent * np.pi / 180) + 0.2)

    [ax.add_patch(w) for w in wedges]
    ax.add_patch(patches.Wedge((0.5, 0), 0.5, 0, 180, color='w'))
    ax.add_patch(patches.Arc((0.5, 0), 1.0, 1.0, theta1=90 - extent, theta2=90 + extent))
    ax.axis('off')
    
    return ax
    