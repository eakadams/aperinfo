#Utility functions

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Handy utility functions
"""

from astropy.io import ascii


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
    apertif_fields = fields[(fields['label'] == 'm') | (fields['label'] == 's') | (fields['label'] == 'l')]
    return apertif_fields['ra'],apertif_fields['dec']
