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
import tol_colors as tc

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
aperinfodir = this_dir[:-9]
filedir = os.path.join(aperinfodir,"files")
figdir = os.path.join(aperinfodir,"figures")

#set up colors - use Paul Tol's  color blind
plt.rc('axes', prop_cycle=plt.cycler('color', list(tc.tol_cset('bright'))))
plt.cm.register_cmap('YlOrBr', tc.tol_cmap('YlOrBr'))
plt.rc('image', cmap='YlOrBr')
prop_cycle = plt.rcParams['axes.prop_cycle']
mpcolors = prop_cycle.by_key()['color']

def plot_sky_view(ra_array_lists, dec_array_lists,
                  label_list, viewname,
                  surveypointings = None,
                  alphalist = None, colorlist = None,
                  show_mds = False, mds_color = None,
                  reobs_ra = None, reobs_dec = None,
                  reobs_color = None):
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
    alphalist : list-like
        List of alpha values (floats) to use in plotting
    colorlist : list-lik
        Optional list of colors to use in plotting
    show_mds : Booelean
        Toggle for whether to outline MDS regions
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
  

    #figure instance
    fig = plt.figure(figsize=(12,12))
    frame = fig.add_axes((0.1,0.1,0.8,0.8))
    f = maputils.FITSimage(externalheader=header)

    #zoom-in to relevant part of graph
    f.set_limits(pxlim = (5,35), pylim=(5,35))
    
    X = np.arange(0,360.0,15.0)
    Y = np.arange(15,90,15) #[20, 30,45, 60, 75]

    
    lon_world = np.arange(0,360,30)
    lat_world = np.arange(0,90,15) #[20, 30, 60, 90]
    lon_constval = 9 * 15 #plo dec labels at 9hours
    lat_constval = 22  #plot ra labels at 22deg

    #set limits to zoom in on relevant areas

    annim = f.Annotatedimage(frame)
    grat = annim.Graticule(axnum=(1,2),wylim=(15.0,90.0), wxlim=(0,360),
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
    ms = 10

    
    #add data
    #check for optional lists
    if alphalist is None:
        alphalist = np.full(len(ra_array_lists),1.0)
    if colorlist is None:
        colorlist = mpcolors[ 0 : len(ra_array_lists)]
    for (ra,dec,lab, alph, color) in zip(ra_array_lists, dec_array_lists,
                                         label_list, alphalist, colorlist):
        xp,yp = annim.topixel(ra,dec)
        annim.Marker(x=xp,y=yp,
                     marker='o',mode='pixel',markersize=ms,
                     color = color, label = lab,
                     alpha = alph)

    #add survey pointings (if provided)
    if surveypointings is not None:
        sra, sdec = get_survey_ra_dec(surveypointings)
        xs,ys = annim.topixel(sra,sdec)
        annim.Marker(x=xs,y=ys,
                     marker='o',mode='pixel',markersize=ms+1,
                     color='black',fillstyle='none')

    #add mds footprints, if specified
    if show_mds is True:
        mds_points = os.path.join(filedir,'mds_pointings.txt')
        mra, mdec = get_survey_ra_dec(mds_points)
        xm, ym = annim.topixel(mra, mdec)
        if mds_color is None:
            mds_color = mpcolors[len(ra_array_lists)]
        annim.Marker(x=xm,y=ym,
                     marker='o',mode='pixel',markersize=ms+1,
                     color=mds_color,fillstyle='none',
                     label = 'Medium-deep',
                     markeredgewidth=2.5)

    #add fields to reobserve, if specified
    if (reobs_ra is not None) and (reobs_dec is not None):
        xr, yr = annim.topixel(reobs_ra, reobs_dec)
        if reobs_color is None:
            reobs_color = mpcolors[len(ra_array_lists)+1]
        annim.Marker(x = xr, y = yr,
                     marker = 'o', mode='pixel', markersize = ms+1,
                     colord = reobs_color,fillstyle = 'none',
                     label = 'To reobserve', markeredgewidth = 2.5)





    #make figure
    annim.plot()
    #add legend
    plt.legend()
    #save fig
    figname = os.path.join(figdir,"{}.pdf".format(viewname))
    plt.savefig(figname)
    plt.close()

     
