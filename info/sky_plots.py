#Functions for plotting on sky

from __future__ import print_function

__author__="E.A.K. Adams"

"""
Functionality for making sky 
overview plots of data
Base files are in "files" directory

Try a class-based approach where I
create a class (centered around an
astropy Table) for each type of data/
table I'm interested in plotting/working
with
"""

import glob
import sys
from astropy.table import Table
import numpy as np
from astropy.io import ascii
from datetime import datetime
import os
import shutil
import matplotlib.pyplot as plt

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
aperinfodir = this_dir[:-4]
filedir = os.path.join(aperinfodir,"files")
#print(filedir)


class ValidCat(object):
    def __init__(self, validfile = os.path.join(filedir,'combined_valid.csv'),
                 obsfile = os.path.join('obs_atdb.csv')):
        """
        """
