#Functions for plotting

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Supporting functionality for making plots
"""

import matplotlib.pyplot as plt
from kapteyn import maputils
#from kapteyn.wcs import galactic, equatorial, fk4_no_e, fk5
from functions.utilities import get_survey_ra_dec
import os
import numpy as np

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
aperinfodir = this_dir[:-9]
filedir = os.path.join(aperinfodir,"files")
figdir = os.path.join(aperinfodir,"figures")

#get mpl colors
prop_cycle = plt.rcParams['axes.prop_cycle']
mpcolors = prop_cycle.by_key()['color']

def plot_sky_view(ra_array_lists, dec_array_lists,
                  label_list, viewname,
                  surveypointings = None):
    """
    Make a sky view plot

    Parameters
    ----------
    ra_array_list : list of arrays
        List of array-like RA values
    dec_array_list : list of arrays
        List of array-like object of Dec values
    label_list : list
        List of string labels
    viewname : str
        Name of view, used in labeling/output
    """

    #start the figure
    #first define a header
    dec0 = 89.9999999999   # Avoid plotting on the wrong side
    header = {'NAXIS'  :  2,
              'NAXIS1' :  40, 'NAXIS2': 40,
              'CTYPE1' : 'RA---ARC',
              'CRVAL1' :  0.0, 'CRPIX1' : 20, 'CUNIT1' : 'deg',
              'CDELT1' : -5.0, 'CTYPE2' : 'DEC--ARC',
              'CRVAL2' :  dec0, 'CRPIX2' : 20,
              'CUNIT2' : 'deg', 'CDELT2' : 5.0,
              }
    """    dec0 = 89.9999999999   # Avoid plotting on the wrong side
    header = {'NAXIS'  :  2,
              'NAXIS1' :  100, 'NAXIS2': 80,
              'CTYPE1' : 'RA---TAN',
              'CRVAL1' :  0.0, 'CRPIX1' : 50, 'CUNIT1' : 'deg',
              'CDELT1' : -5.0, 'CTYPE2' : 'DEC--TAN',
              'CRVAL2' :  dec0, 'CRPIX2' : 40,
              'CUNIT2' : 'deg', 'CDELT2' : 5.0,
              }
    header = {'NAXIS'  :  2,
              'NAXIS1' :  100, 'NAXIS2': 80,
              'CTYPE1' : 'RA---ARC',
              'CRVAL1' :  0.0, 'CRPIX1' : 50, 'CUNIT1' : 'deg',
              'CDELT1' : -5.0, 'CTYPE2' : 'DEC--ARC',
              'CRVAL2' :  dec0, 'CRPIX2' : 40,
              'CUNIT2' : 'deg', 'CDELT2' : 5.0,
              }
    """
  
    X = np.arange(0,360.0,15.0)
    Y = np.arange(0,90,15) #[20, 30,45, 60, 75]


    #figure instance
    fig = plt.figure(figsize=(12,12))
    frame = fig.add_axes((0.1,0.1,0.8,0.8))
    f = maputils.FITSimage(externalheader=header)

    lon_world = np.arange(0,360,30)
    lat_world = np.arange(0,90,15) #[20, 30, 60, 90]
    lon_constval = None
    lat_constval = 20

    annim = f.Annotatedimage(frame)
    grat = annim.Graticule(axnum=(1,2),wylim=(0.0,90.0), wxlim=(0,360),
                       startx=X, starty=Y)

    grat.setp_gratline(color='0.75')
    
    grat.setp_lineswcs1(20, color='g', linestyle='--')

    # Plot labels inside the plot
    lon_kwargs = {'color':'k', 'fontsize':12}
    lat_kwargs = {'color':'k', 'fontsize':12}
    grat.Insidelabels(wcsaxis=0,
                      world=lon_world, constval=lat_constval,
                      fmt="Hms",
                      **lon_kwargs)
    grat.Insidelabels(wcsaxis=1,
                      world=lat_world, constval=lon_constval,
                      fmt="Dms",
                      **lat_kwargs)

    #set marker size
    ms = 8

    #add survey pointings (if provided)
    if surveypointings is not None:
        sra, sdec = get_survey_ra_dec(surveypointings)
        xs,ys = annim.topixel(sra,sdec)
        annim.Marker(x=xs,y=ys,
                     marker='o',mode='pixel',markersize=ms,
                     color='black',fillstyle='none')

    #add data
    for i,(ra,dec,lab) in enumerate(zip(ra_array_lists, dec_array_lists,
                                        label_list)):
        xp,yp = annim.topixel(ra,dec)
        annim.Marker(x=xp,y=yp,
                     marker='o',mode='pixel',markersize=ms,
                     color = mpcolors[i], label = lab)

    #make figure
    annim.plot()
    #add legend
    plt.legend()
    #save fig
    figname = os.path.join(figdir,"census_{}.pdf".format(viewname))
    plt.savefig(figname)
    plt.close

     
