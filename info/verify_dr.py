#Verify data release
from __future__ import print_function

__author__="E.A.K. Adams"

"""
Functions to check data release tables
pulled from VO services against provided
files
"""

import glob
import sys
from astropy.table import Table
import numpy as np
from astropy.io import ascii
import os
import pyvo
import astropy.io.fits as fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import requests


#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
aperinfodir = this_dir[:-4]
filedir = os.path.join(aperinfodir,"files")
tabledir = os.path.join(aperinfodir,"tables")
#print(filedir)


def get_number_ind_beams():
    """
    Get the number of independent beams released
    E.g., amount of unique sky coverage
    """
    dr_csv = ascii.read(os.path.join(tabledir,"dr_year1_cont.csv"))
    print(dr_csv.colnames)
    names_dr = np.array( ["{0}-{1}".format(field,bm) for field,bm in dr_csv['Name','Beam']] )
    print(len(names_dr),len(np.unique(names_dr)))
    
def check_vo_site_cont():
    """
    Check the continuum table
    retrieved as csv directly from 
    vo.astron.nl
    """
    vo_csv = ascii.read(os.path.join(filedir,"table_cont_vo.csv"))
    dr_csv = ascii.read(os.path.join(tabledir,"dr_year1_cont.csv"))
    #print(vo_csv.colnames)
    #print(dr_csv.colnames)
    #Get ObsID / beam names and match those (should be unique)
    names_vo = np.array( ["{0}-{1}".format(obsid,bm) for obsid,bm in vo_csv['obsid','beam_number']] )
    #print(names_vo[0:10])
    names_dr = np.array( ["{0}-{1}".format(obsid,bm) for obsid,bm in dr_csv['ObsID','Beam']] )
    #print(names_dr[0:10])
    #check for names in one but not the other
    missing_vo = np.setdiff1d(names_dr,names_vo)
    extra_vo = np.setdiff1d(names_vo,names_dr)
    print(missing_vo)
    print(extra_vo)
    print(len(dr_csv))


def check_vo_auto():
    """
    Check what is in VO by retrieving with pyvo
    Also checking pyvo code
    """
    #connect to TAP query
    tap_service = pyvo.dal.TAPService('https://vo.astron.nl/__system__/tap/run/tap')
    for table in tap_service.tables:
        print(table.name)
    cont_table = tap_service.tables['apertif_dr1.continuum_images']
    print(cont_table.columns)
    print('Performing TAP query')
    result = tap_service.search(
        "SELECT TOP 5 target, beam_number, accref, centeralpha, centerdelta, obsid, DISTANCE(" \
        "POINT('ICRS', centeralpha, centerdelta),"\
        "POINT('ICRS', 208.36, 52.36)) AS dist"\
        " FROM apertif_dr1.continuum_images"  \
        " WHERE 1=CONTAINS("
        "    POINT('ICRS', centeralpha, centerdelta),"\
        "    CIRCLE('ICRS', 208.36, 52.36, 0.08333333)) "\
        " ORDER BY dist ASC"
    )
    print(result)
    # The result can also be obtained as an astropy table
    astropy_table = result.to_table()
    print('-' * 10 + '\n' * 3)
    ## You can also download and plot the image
    import astropy.io.fits as fits
    from astropy.wcs import WCS
    import matplotlib.pyplot as plt
    import requests, os
    import numpy as np

    # DOWNLOAD only the first result
    #
    print('Downloading only the first result')
    file_name = '{}_{}_{}.fits'.format(
        result[0]['obsid'].decode(),
        result[0]['target'].decode(),
        result[0]['beam_number'])
    path = os.path.join(os.getcwd(), file_name)
    http_result = requests.get(result[0]['accref'].decode())
    print('Downloading file in', path)
    with open(file_name, 'wb') as fout:
        for content in http_result.iter_content():
            fout.write(content)
    hdu = fits.open(file_name)[0]
    wcs = WCS(hdu.header)
    # dropping unnecessary axes
    wcs = wcs.dropaxis(2).dropaxis(2)
    plt.subplot(projection=wcs)
    plt.imshow(hdu.data[0, 0, :, :], vmax=0.0005)
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.show()


def get_daily_image():
    """
    Use VO tools to get and plot 
    continuum image for daily image
    """
     #connect to TAP query
    tap_service = pyvo.dal.TAPService('https://vo.astron.nl/__system__/tap/run/tap')
    cont_table = tap_service.tables['apertif_dr1.continuum_images']
    print('Performing TAP query')
    result = tap_service.search(
        "SELECT TOP 5 target, beam_number, accref, centeralpha, centerdelta, obsid, DISTANCE(" \
        "POINT('ICRS', centeralpha, centerdelta),"\
        "POINT('ICRS', 239.37, 54.75)) AS dist"\
        " FROM apertif_dr1.continuum_images"  \
        " WHERE 1=CONTAINS("
        "    POINT('ICRS', centeralpha, centerdelta),"\
        "    CIRCLE('ICRS', 239.37, 54.75, 0.08333333)) "\
        " ORDER BY dist ASC"
    )
    print(result)
    # The result can also be obtained as an astropy table
    astropy_table = result.to_table()
    print('-' * 10 + '\n' * 3)
    ## You can also download and plot the image
    

    # DOWNLOAD only the first result
    #
    print('Downloading only the first result')
    file_name = '{}_{}_{}.fits'.format(
        result[0]['obsid'].decode(),
        result[0]['target'].decode(),
        result[0]['beam_number'])
    path = os.path.join(os.getcwd(), file_name)
    http_result = requests.get(result[0]['accref'].decode())
    print('Downloading file in', path)
    with open(file_name, 'wb') as fout:
        for content in http_result.iter_content():
            fout.write(content)
    hdu = fits.open(file_name)[0]
    wcs = WCS(hdu.header)
    # dropping unnecessary axes
    wcs = wcs.dropaxis(2).dropaxis(2)
    plt.subplot(projection=wcs)
    plt.imshow(hdu.data[0, 0, :, :], vmax=0.0005)
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.show()

def plot_daily_image():
    """
    Having retrieved file, fuss with plotting here
    """
    #connect to TAP query
    tap_service = pyvo.dal.TAPService('https://vo.astron.nl/__system__/tap/run/tap')
    cont_table = tap_service.tables['apertif_dr1.continuum_images']
    print('Performing TAP query')
    result = tap_service.search(
        "SELECT TOP 5 target, beam_number, accref, centeralpha, centerdelta, obsid, DISTANCE(" \
        "POINT('ICRS', centeralpha, centerdelta),"\
        "POINT('ICRS', 239.37, 54.75)) AS dist"\
        " FROM apertif_dr1.continuum_images"  \
        " WHERE 1=CONTAINS("
        "    POINT('ICRS', centeralpha, centerdelta),"\
        "    CIRCLE('ICRS', 239.37, 54.75, 0.08333333)) "\
        " ORDER BY dist ASC"
    )
    print(result)
    # The result can also be obtained as an astropy table
    astropy_table = result.to_table()
    print('-' * 10 + '\n' * 3)
    ## You can also download and plot the image
    

    # DOWNLOAD only the first result
    #
    print('Downloading only the first result')
    file_name = '{}_{}_{}.fits'.format(
        result[0]['obsid'].decode(),
        result[0]['target'].decode(),
        result[0]['beam_number'])
    hdu = fits.open(file_name)[0]
    wcs = WCS(hdu.header)
    # dropping unnecessary axes
    wcs = wcs.dropaxis(2).dropaxis(2)
    print(hdu.data.shape)
    plt.subplot(projection=wcs)
    #plt.subplot(projection=wcs)
    plt.imshow(hdu.data[0, 0, 1000:2000, 1000:2000], vmax=0.0005)
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.show()

