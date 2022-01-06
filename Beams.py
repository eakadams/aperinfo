#Beams object

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Object to hold Apertif beam info
Children objects for more specific purposes:
 - DR1 : First data release paper
"""

import os
import matplotlib.pyplot as plt
from astropy.io import ascii
import numpy as np
from astropy.table import Table, join
from functions.plots import plot_sky_view

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
#at a differnt level this time, so actually in this dir
aperinfodir = this_dir
filedir = os.path.join(aperinfodir,"files")
tabledir = os.path.join(aperinfodir,"tables")
figdir = os.path.join(aperinfodir,"figures")


#get mpl colors
prop_cycle = plt.rcParams['axes.prop_cycle']
mpcolors = prop_cycle.by_key()['color']


#set fonts and latex
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})


class Beams(object):
    """
    A class holding Apertif beam-based information

    ...

    Attributes
    ----------
    obsinfo : Astropy Table
         table of observational information
    beaminfo : Astropy Table
         table of beam-based information

    Methods
    -------
    
    """
    def __init__(self,
                 obsfile = os.path.join(filedir,'obsatdb.csv'),
                 happilifile = os.path.join(filedir,'happili.csv'),
                 contfile = os.path.join(filedir,"cont_allbeams.csv"),
                 hifile = os.path.join(filedir,"hi_allbeams.csv"),
                 polfile = os.path.join(filedir,"pol_allbeams.csv")
                 ):
        """
        Construct initial attributes for Beam object

        This includes combining the validation information for three
        separate data products: continuum, polarization, line (HI)

        Parameters:
        -----------
        obsfile : str
            Full path to atdb observational file
        happilifile : str
            Full path to file with happili processing information
        contfile : str
            Full path to continuum validation file
        hifile : str
            Full path to HI validation file
        polfile : str
            Full path to polarization validation file
        """

        self.obsinfo = ascii.read(obsfile)

        #setup beam info based on obsinfo
        beam_table_length = len(self.obsinfo) * 40
        self.beaminfo = Table()
        self.beaminfo['ObsID'] = np.empty(beam_table_length, dtype = int)
        self.beaminfo['beam'] = np.empty(beam_table_length, dtype = int)
        for n,id in enumerate(self.obsinfo['taskID']):
            ind = n*40
            self.beaminfo['ObsID'][ind:ind+40] = np.full(40,id)
            self.beaminfo['beam'][ind:ind+40] = np.arange(40)

            
