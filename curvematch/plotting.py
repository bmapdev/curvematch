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


def plotcurve(q1, interactive=False, filename=None):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    if q1.dim == 3:
        ax.plot(q1.coords[0, :], q1.coords[1, :], q1.coords[2, :])
    else:
        ax.plot(q1.coords[0, :], q1.coords[1, :])

    ax.view_init(90, -90)
    if interactive:
        plt.show()
    else:
        if filename is None:
            filename = 'curve_plot'
        plt.savefig(path.join(Settings.output_dir, filename))


def path_plot():
    pass


def plot_deformationfield():
    pass

