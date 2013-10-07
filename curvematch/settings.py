#!/usr/local/epd/bin/python

"""Settings and Configurations"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"


class Settings(object):

    def __init__(self, interactive_mode=False, backend='PDF'):
        self.interactive_mode = interactive_mode
        self.backend = backend

    interactive_mode = False
    backend = 'PDF'
    output_dir = '.'
