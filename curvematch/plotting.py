#!/usr/local/epd/bin/python
"""Plotting functions"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"


from settings import Settings
import matplotlib
if not Settings.interactive_mode:
    matplotlib.use(Settings.backend)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from os import path


def plot_curve(q1, fig_num=1, interactive=False, filename=None):
    fig = plt.figure(fig_num)
    ax = fig.gca(projection='3d')
    if q1.dim == 3:
        ax.plot(q1.coords[0, :], q1.coords[1, :], q1.coords[2, :])
    else:
        ax.plot(q1.coords[0, :], q1.coords[1, :])

    plt.autoscale(tight=True)
    plt.axis('equal')
    plt.axis('off')
    ax.view_init(90, -90)

    # plt.tight_layout(0.0)
    # fig.set_tight_layout(True)
    # plt.tight_layout()
    if interactive:
        plt.show()
    else:
        if filename is None:
            filename = 'curve_plot'
        # SaveFigureAsImage(path.join(Settings.output_dir, filename), fig)
        plt.savefig(path.join(Settings.output_dir, filename), bbox_inches='tight', pad_inches=0.0, dpi=300)
    return fig

def plot_path(curve_path, filename=None):
    for i in xrange(0, len(curve_path)):
        curve_path[i].coords[0,:] += 0.2*i
        fig = plot_curve(curve_path[i], interactive=False, filename=filename)
        fig.hold()


def plot_deformationfield():
    pass

def SaveFigureAsImage(fileName,fig=None,**kwargs):
    ''' Save a Matplotlib figure as an image without borders or frames.
       Args:
            fileName (str): String that ends in .png etc.

            fig (Matplotlib figure instance): figure you want to save as the image
        Keyword Args:
            orig_size (tuple): width, height of the original image used to maintain
            aspect ratio.
    '''
    fig_size = fig.get_size_inches()
    w,h = fig_size[0], fig_size[1]
    fig.patch.set_alpha(0)
    if kwargs.has_key('orig_size'): # Aspect ratio scaling if required
        w,h = kwargs['orig_size']
        w2,h2 = fig_size[0],fig_size[1]
        fig.set_size_inches([(w2/w)*w,(w2/w)*h])
        fig.set_dpi((w2/w)*fig.get_dpi())
    a=fig.gca()
    a.set_frame_on(False)
    a.set_xticks([]); a.set_yticks([])
    plt.axis('off')
    plt.xlim(0,h); plt.ylim(w,0)
    fig.savefig(fileName, transparent=True, bbox_inches='tight', \
                        pad_inches=0)