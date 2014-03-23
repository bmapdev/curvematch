#!/usr/local/epd/bin/python
"""Test Curve Matching functions"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from curvematch.settings import Settings
from shapeio import curveio
from curvematch.qshape import QShape
import geodesics
from curve import Curve
import numpy as np
from math import pi
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import os



def match_curve_pair(curvefilename1, curvefilename2, settings, rotation=True, siz=100, return_curves=False):

    c1 = Curve(file=curvefilename1)
    c2 = Curve(file=curvefilename2)
    c1.resample_curve_uniform(siz)
    c2.resample_curve_uniform(siz)

    q1 = QShape()
    q2 = QShape()
    q1.from_curve(c1)
    q2.from_curve(c2)

    if settings.closed:
        geodesic = geodesics.compute_for_closed_curves(q1, q2, settings)
    else:
        geodesic = geodesics.compute_for_open_curves_elastic(q1, q2, settings, rotation)

    if return_curves:
        return geodesic, c1, c2
    else:
        return geodesic


def elastic_curve_matching(template_curve, match_curve, settings, rotation=True):
    """
    Finds finds matching between curves using dynamic programing and,
    applies the gamma values from this matching to elastically match
    match_curve to template, and returns the matched curve.
    """
    qt = QShape()
    qm = QShape()
    qt.from_curve(template_curve)
    qm.from_curve(match_curve)
    if settings.closed:
        geodesic = geodesics.compute_for_closed_curves_elastic(qt, qm, settings)
    else:
        geodesic = geodesics.compute_for_open_curves_elastic(qt, qm, settings, rotation)

    matched_curve = match_curve
    for i in xrange(match_curve.dim()):
        matched_curve.coords[i, :] = np.interp(geodesic.gamma, np.linspace(0, 2*pi, match_curve.siz()),
                                               match_curve.coords[i, :])
    match_curve.gamma = geodesic.gamma
    return matched_curve


def get_group_matching(template_curve, match_curves, settings=False, resample_size=500, rotation=True, match=True):
    """
    Resamples all curves, and then preforms an elastic matching
    between every curve in match_curve_list and the template_curve.
    Returns a list of adjusted curves with gamma values as attributes.
    """
    if not settings:
        settings = geodesics.Geodesic()
        settings.steps = 5
        settings.closed = True

    if type(match_curves) != list and type(match_curves) != tuple:
        match_curves = [match_curves]
    template_curve.resample_curve_uniform(resample_size)
    matched_curves_list = []
    for current_curve in match_curves:
        current_curve.resample_curve_uniform(resample_size)
        if match:
            matched_curves_list.append(elastic_curve_matching(template_curve, current_curve, settings, rotation))
        else:
            matched_curves_list.append(current_curve)
    return matched_curves_list


def load_curves_from_top_and_bottom(top_curves_file, bot_curves_file, connect=False):
    top_paths = open(top_curves_file, 'r')
    bot_paths = open(bot_curves_file, 'r')
    top_paths = top_paths.read().split('\n')
    bot_paths = bot_paths.read().split('\n')

    if len(top_paths) != len(bot_paths):
        raise ValueError("Unequal number of top and bottom curves!")

    top_curves = []
    bot_curves = []
    connected_curves = []
    for i in xrange(len(top_paths)):
        print '\n', top_paths[i], '\n'
        current_top = Curve(file=top_paths[i])
        current_bot = Curve(file=bot_paths[i])
        if not connect:
            top_curves.append(current_top)
            current_bot.append(current_bot)
        else:
            #print current_top.coords ,'\n'
            current_top.append_curve(current_bot)
            connected_curves.append(current_top)
            #print current_top.coords ,'\n'
    if not connect:
        return [top_curves, bot_curves]
    else:
        return connected_curves


def group_matching_batch(top_curves_file, bot_curves_file):
    """
    Reads in paths from newline seperated text files. The first curve listed will be
    used as the template curve.
    """
    curves_list = load_curves_from_top_and_bottom(top_curves_file, bot_curves_file, connect=True)
    template = curves_list[0]
    #print template
    settings = geodesics.Geodesic()
    settings.steps = 5
    settings.closed = True

    i = 1
    matched_curves = get_group_matching(template, curves_list[1:], match=False)
    for curve in matched_curves:
        plot_matching("subject_uniform"+str(i), template, curve)
        i += 1

    matched_curves = get_group_matching(template, curves_list[1:], match=True)
    i = 1
    for curve in matched_curves:
        plot_matching("subject_closed"+str(i), template, curve)
        i += 1
    settings.closed = False

    i = 1
    for curve in matched_curves:
        plot_matching("subject_open"+str(i), template, curve)
        i += 1



def plot_matching(plot_title, curve1, curve2, lines=20, offset=5):

    dims = [0,1,2]
    dims.remove(curve2.least_variant_dimension())

    font = 11
    shift_x = np.min(curve1.coords[dims[0], :]) - np.min(curve2.coords[dims[0], :])
    shift_y = np.max(curve1.coords[dims[1], :]) - np.max(curve2.coords[dims[1], :])

    curve2.coords[dims[0]] += shift_x
    curve2.coords[dims[1]] += shift_y
    shift_y_down = max(curve2.coords[dims[0]]) - min(curve1.coords[dims[1]]) + offset
    curve2.coords[dims[1]] -= shift_y_down
    plt.tick_params(labelbottom=False, labelleft=False, labelright=False)
    plt.subplot(111)
    if plot_title:
        plt.title(plot_title, fontsize=font)
    line_step = curve1.siz() / lines
    plt.plot(curve1.coords[dims[0]], curve1.coords[dims[1]])
    plt.plot(curve2.coords[dims[0]], curve2.coords[dims[1]])
    for i in xrange(0, curve1.siz(), line_step):
        plt.plot([curve1.coords[dims[0]][i], curve2.coords[dims[0]][i]],
                 [curve1.coords[dims[1]][i], curve2.coords[dims[1]][i]])
    plt.savefig(plot_title + '.pdf') ##Just showing for testing purposes
    print "Plot for Subject ", plot_title," has been saved to ", os.getcwd()
    plt.close('all')




































