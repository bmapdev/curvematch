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
import plotting
import os


def match_curve_pair(curve_target_filename, curve_source_filename, settings, rotation=True, siz=100, return_curves=False, linear=False):

    c_target = Curve(file=curve_target_filename)
    c_source = Curve(file=curve_source_filename)
    c_target.resample_curve_uniform(siz)
    c_source.resample_curve_uniform(siz)

    q_target = QShape()
    q_source = QShape()
    q_target.from_curve(c_target)
    q_source.from_curve(c_source)

    if settings.closed:
        geodesic = geodesics.compute_for_closed_curves(q_target, q_source, settings)
    else:
        geodesic = geodesics.compute_for_open_curves_elastic(q_target, q_source, settings, rotation, linear)

    if return_curves:
        c_source_matched = c_source.return_reparameterized_by_gamma(geodesic.gamma)
        return geodesic, c_target, c_source_matched
    else:
        return geodesic


def match_curve_group(src_curve_paths, target, openflag=False, linearflag=False, norotateflag=False, resize=200):
    settings = Settings()

    if openflag:
        settings.closed = False

    geodesic_array = []
    src_curve_matched_to_target_array = []
    target_curve_array = []

    for idx, src_curve in enumerate(src_curve_paths):
        geodesic, target_curve, src_curve_matched_to_target = match_curve_pair(target, src_curve, settings, rotation=norotateflag, siz=resize, return_curves=True, linear=linearflag)
        geodesic_array.append(geodesic)
        src_curve_matched_to_target_array.append(src_curve_matched_to_target)
        target_curve_array.append(target_curve)

    return geodesic_array, src_curve_matched_to_target_array, target_curve_array


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
    match_curve.geodesic = geodesic
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
    #standardize_callosal_curve_coordinates(template_curve)
    template_curve.resample_curve_uniform(resample_size)
    matched_curves_list = []
    for current_curve in match_curves:
        #standardize_callosal_curve_coordinates(current_curve)
        current_curve.resample_curve_uniform(resample_size)
        if match:
            matched_curves_list.append(elastic_curve_matching(template_curve, current_curve, settings, rotation))
        else:
            matched_curves_list.append(current_curve)
    return matched_curves_list


def load_curves_from_top_and_bottom(top_curves_file, bot_curves_file, connect=False):
    """
    Reads in paths from either from
    1) Newline seperated text files.
        e.g. group_match_batch("~/top_curve_paths.txt,~/bot_curve_paths.txt")
    2) Ordered lists of curve file paths.
        e.g. group_match_batch([top_curve_path_1,top_curve_path_2,...],
                               [bot_curve_path_1,bot_curve_path_2,...])

     The first curve listed will be
    used as the template curve.
    """
    if type(top_curves_file) == str:
        top_paths = open(top_curves_file, 'r')
        top_paths = top_paths.read().split('\n')
    else:
        top_paths = top_curves_file

    if type(bot_curves_file) == str:
        bot_paths = open(bot_curves_file, 'r')
        bot_paths = bot_paths.read().split('\n')
    else:
        bot_paths = bot_curves_file

    if len(top_paths) != len(bot_paths):
        raise ValueError("Unequal number of top and bottom curves!")

    top_curves = []
    bot_curves = []
    connected_curves = []
    for i in xrange(len(top_paths)):
        print '\n', top_paths[i], '\n'
        current_top = Curve(file=top_paths[i])
        current_bot = Curve(file=bot_paths[i])
        standardize_callosal_curve_coordinates(current_top)
        standardize_callosal_curve_coordinates(current_bot)
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

    curves_list = load_curves_from_top_and_bottom(top_curves_file, bot_curves_file, connect=True)
    template = curves_list[0]
    #print template
    settings = geodesics.Geodesic()
    settings.steps = 5
    settings.closed = True

    i = 1
    matched_curves = get_group_matching(template, curves_list[1:], match=False)
    for curve in matched_curves:
        plotting.plot_matching("uniform_match"+str(i), template, curve)
        i += 1

    matched_curves = get_group_matching(template, curves_list[1:], match=True)
    i = 1
    for curve in matched_curves:
        plotting.plot_matching("elastic_closed_match"+str(i), template, curve)
        i += 1
    settings.closed = False

    i = 1
    for curve in matched_curves:
        plotting.plot_matching("elastic_open_match"+str(i), template, curve)
        i += 1
    i = 1
    for curve in matched_curves:
        plotting.simple_curve_plot("gamma"+str(i), curve.gamma)
        i += 1


def standardize_callosal_curve_coordinates(curve):
    varance_ditionary = {}
    for dim in xrange(curve.dim()):
        varance_ditionary[np.std(curve.coords[dim, :])] = dim
    keys = sorted(varance_ditionary.keys(), reverse=True)

    newCoords = np.zeros((curve.dim(), curve.siz()))
    i = 0
    for key in keys:
        newCoords[i, :] = curve.coords[varance_ditionary[key], :]
        i += 1
    #x = 0
   # print curve.coords[x,0], " -> ",curve.coords[x,-1],'\n'
    #if newCoords[x,0] > newCoords[x, -1]:
     #   print "YUP!\n"
      #  newCoords = newCoords[:, ::-1]

    curve.coords = newCoords
