__author__ = "Brandon R. Ayers"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from curvematch.test.test_DP_matchcy import test_DP_matchcy
from curvematch.settings import Settings
from curvematch.qshape import QShape
from shapeio import curveio
import pkg_resources
import sys


def run_test(test_id):
    if test_id == 1:
        return test_DP_matchcy()

    elif test_id == 2:
        pass
if __name__ == '__main__':
    run_test(sys.argv[1])



