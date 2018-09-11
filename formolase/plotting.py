def box_off(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

def deduplicate_legend_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
        if label not in newLabels:
            newLabels.append(label)
            newHandles.append(handle)
    return (newHandles, newLabels)

def set_matplotlib_font_sizes(
    small_size = 12,
    medium_size = 14,
    large_size = 16,):
    import matplotlib.pyplot as plt

    plt.rc('font', size=small_size)          # controls default text sizes
    plt.rc('axes', titlesize=large_size)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium_size)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small_size)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small_size)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small_size)    # legend fontsize
    plt.rc('figure', titlesize=large_size)  # fontsize of the figure title

def set_jue_plotting_defaults():
    import matplotlib.pyplot as plt
    set_matplotlib_font_sizes()
    plt.style.use('default')
    plt.style.use('default')
    
def set_outside_axes_labels(ax, axs, xlabel, ylabel):
    if any([ax==axs[i,0] for i in range(axs.shape[0])]):
        ax.set_ylabel(ylabel);
    if any([ax==axs[axs.shape[0]-1,i] for i in range(axs.shape[1])]):
        ax.set_xlabel(xlabel);

def savefig(fig, filename, **kwargs):
    args = {'dpi':150,
            'bbox_inches':'tight'
           }
    args.update(kwargs)
    fig.savefig(filename, **args)
