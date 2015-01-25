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
import os
import numpy as np


def plot_curve(q1, fig_num=1, interactive=False, filename=None):
    fig = plt.figure(fig_num)
    ax = fig.gca(projection='3d')
    if q1.dim() == 3:
        ax.plot(q1.coords[0, :], q1.coords[1, :], q1.coords[2, :])
    else:
        ax.plot(q1.coords[0, :], q1.coords[1, :])

    plt.autoscale(tight=True)
    plt.axis('equal')
    plt.axis('off')
    # ax.view_init(0, 90)
    ax.autoscale_view()

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


def SaveFigureAsImage(fileName, fig=None, **kwargs):
    """ Save a Matplotlib figure as an image without borders or frames.
       Args:
            fileName (str): String that ends in .png etc.

            fig (Matplotlib figure instance): figure you want to save as the image
        Keyword Args:
            orig_size (tuple): width, height of the original image used to maintain
            aspect ratio.
    """


    fig_size = fig.get_size_inches()
    w, h = fig_size[0], fig_size[1]
    fig.patch.set_alpha(0)
    if 'orig_size' in kwargs: # Aspect ratio scaling if required
        w,h = kwargs['orig_size']
        w2,h2 = fig_size[0],fig_size[1]
        fig.set_size_inches([(w2/w)*w,(w2/w)*h])
        fig.set_dpi((w2/w)*fig.get_dpi())
    a = fig.gca()
    a.set_frame_on(False)
    a.set_xticks([])
    a.set_yticks([])
    plt.axis('off')
    plt.xlim(0, h)
    plt.ylim(w, 0)
    fig.savefig(fileName, transparent=True, bbox_inches='tight', \
                pad_inches=0)


def plot_matching(plot_title, curve1, curve2, lines=30, offset=5, outdir=''):

    dims1 = [0, 1, 2]
    dims1.remove(curve1.least_variant_dimension())

    dims2 = [0, 1, 2]
    dims2.remove(curve2.least_variant_dimension())
    font = 11
    shift_x = np.min(curve1.coords[dims1[0], :]) - np.min(curve2.coords[dims2[0], :])
    shift_y = np.max(curve1.coords[dims1[1], :]) - np.max(curve2.coords[dims2[1], :])

    curve2.coords[dims2[0]] += shift_x
    curve2.coords[dims2[1]] += shift_y
    shift_y_down = max(curve2.coords[dims2[0]]) - min(curve1.coords[dims1[1]]) + offset
    curve2.coords[dims2[1]] -= shift_y_down
    plt.tick_params(labelbottom=False, labelleft=False, labelright=False)
    plt.subplot(111)
    if plot_title:
        plt.title(plot_title, fontsize=font)
    line_step = curve1.siz() / lines
    plt.plot(curve1.coords[dims1[0]], curve1.coords[dims1[1]])
    plt.plot(curve2.coords[dims2[0]], curve2.coords[dims2[1]])
    for i in xrange(0, curve1.siz(), line_step):
        plt.plot([curve1.coords[dims1[0]][i], curve2.coords[dims2[0]][i]],
                 [curve1.coords[dims1[1]][i], curve2.coords[dims2[1]][i]])
    plt.savefig(os.path.join(outdir, plot_title) + '.pdf')
    print "Plot for Subject ", plot_title, " has been saved"
    plt.close('all')


def scalar_function_plot(plot_title, coords):
    plt.tick_params(labelbottom=False, labelleft=False, labelright=False)
    plt.subplot(111)
    if plot_title:
        plt.title(plot_title, fontsize=12)
    plt.plot(coords)
    plt.savefig(plot_title + '.pdf')
    print "Plot for Subject ", plot_title ," has been saved"
    plt.close('all')
