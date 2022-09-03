#Utility functions

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Handy utility functions
"""

from astropy.io import ascii
from astropy.coordinates import SkyCoord
import numpy as np
import os
import astropy.units as u

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
aperinfodir = this_dir[:-9]
filedir = os.path.join(aperinfodir,"files")
#print(filedir)

def get_med_onesig(values):
    """
    Get median and one-sigma percentiles

    Parameters:
    - values : nparray
    
    Returns:
    median : float
    per16 : 16th percentile
    per84 : 84th percentile
    """
    median = np.nanmedian(values)
    per16 = np.nanpercentile(values, 16)
    per84 = np.nanpercentile(values, 84)

    lower_onesig = median - per16
    upper_onesig = per84 - median
    
    return median, lower_onesig, upper_onesig
    

def get_beam_ra_dec(field_ra, field_dec, beams,
                    beampos =  os.path.join(filedir,'cb_offsets.txt')):
    """
    For matched arrays of beams and field center ra/decs, get the
    ra, dec center of the specific beam

    Parameters:
    -----------
    - field_ra : array-like, R.A. of the field in degrees
    - field_dec : array-like, Dec. of the field in degrees
    - beams : array-like, compound beam
    - beampos : path to compound beam position file

    Returns:
    --------
    - ra : np-array of RA in degrees
    - dec : np-array of Dec in degrees
    """
    #get field skycoord object
    field_coords = SkyCoord(field_ra, field_dec, frame='icrs',unit='deg')                

    #read in CB positions
    cbpos = ascii.read(beampos)

    #set up ra,dec output
    ra = np.empty(len(field_ra))
    dec = np.empty(len(field_ra))

    #iterate through each entry
    #think there has to a cleaner way to do this, but I'm lazy
    for i, (coord, bm) in enumerate(zip(field_coords, beams)):
        #get cb position for beam
        cb_ra_off = cbpos['ra'][bm]
        cb_dec_off = cbpos['dec'][bm]

        #then add offsets to field position, accounting for declination
        dec_beam = coord.dec.deg + cb_dec_off
        ra_beam = coord.ra.deg + ( cb_ra_off /
                                   np.cos(coord.dec.to(u.radian).value) )

        ra[i] = ra_beam
        dec[i] = dec_beam

    return ra, dec    


def get_survey_ra_dec(survey_pointings):
    """
    Take survey pointing file and return ra and dec array

    Parameters
    ----------
    survey_pointings : str
        Full path to file of survey pointings

    Returns
    -------
    ra : array
        Array of decimal degree RA values
    dec : array
        Array of decimal degree Dec values
    """
    #read  in file
    fields = ascii.read(survey_pointings,format='fixed_width')
    #find only survey pointings
    apertif_fields = fields[(fields['label'] == 'm') | (fields['label'] == 's') |
                            (fields['label'] == 'l')  | (fields['label'] == 'd')]
    return apertif_fields['ra'],apertif_fields['dec']
