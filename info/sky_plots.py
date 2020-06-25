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
from astropy.coordinates import SkyCoord
from astropy import units as u
from kapteyn import maputils
from kapteyn.wcs import galactic, equatorial, fk4_no_e, fk5

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
aperinfodir = this_dir[:-4]
filedir = os.path.join(aperinfodir,"files")
#print(filedir)


"""
Classes for different types of tables
"""

class ValidCat(object):
    def __init__(self, validfile = os.path.join(filedir,'combined_valid.csv'),
                 obsfile = os.path.join(filedir,'obsatdb.csv'),
                 beampos = os.path.join(filedir,'cb_offsets.txt')):
        """
        Initalize ValidCat object
        This has the validation information
        Inputs:
        - validfile(str): File with validation information. Default to standard w/in package
        - obsfile (str): File with observation information (taskid, ra/dec). Default to standard w/in package
        - beampos (str): File with compound beam offsets. Default to standard w/in package
        """
        #master table of ValidCat
        self.table = ascii.read(validfile)
        #obsinfo as supplementary information
        self.obsinfo = ascii.read(obsfile)
        #and beam position
        self.cb_pos = ascii.read(beampos)

        #as part of initialization, I want to get coords of every entry in table
        #use a helper function defined elsewhere for this
        ra_array, dec_array = get_ra_dec(self.table['taskid'],
                                         self.table['beam'],
                                         self.obsinfo['taskID'],
                                         self.obsinfo['field_ra'],
                                         self.obsinfo['field_dec'],
                                         self.cb_pos)
        self.table['RA'] = ra_array
        self.table['Dec'] = dec_array

    def plot_observed(self,figname=None,kapteyn=True):
        """
        Make a plot of observed fields, build from base code
        """
        if figname is None:
            figname = os.path.join(filedir,"observed_fields.pdf")
        ra_list = [self.obsinfo['field_ra']]
        dec_list = [self.obsinfo['field_dec']]
        color_list = ['blue']
        label_list = ['Observed field']
        if kapteyn is True:
            sky_plot_kapteyn(ra_list,dec_list,color_list,label_list,figname)
        else:
            sky_plot(ra_list,dec_list,color_list,label_list,figname)


"""
Helper functions
"""
def sky_plot_kapteyn(ra_array_lists,dec_array_lists,
                     color_list,label_list,figname):
    """
    Make sky plots, using Kapteyn python package
    Inputs:
    ra_array_lists: list of RA arrays
    dec_array_lists: list of Dec array
    color_list: List of colors for plotting
    label_list: List of labels
    figname: output name for figure
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
    annim = f.Annotatedimage(frame)
    grat = annim.Graticule(axnum=(1,2),wylim=(0.0,90.0), wxlim=(0,360),
                       startx=X, starty=Y)
    #grat = annim.Graticule(wylim=(0.0,90.0), wxlim=(0,360),
    #                       startx=X, starty=Y)
    grat.setp_gratline(color='0.75')
    lon_world = np.arange(0,360,30)
    lat_world = np.arange(0,90,15) #[20, 30, 60, 90]
    grat.setp_lineswcs1(20, color='g', linestyle='--')

    # Plot labels inside the plot
    lon_constval = None
    lat_constval = 20
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
    #iterate through arrays and add to plot
    for ra,dec,cname,labname in zip(ra_array_lists,
                                    dec_array_lists,color_list,label_list):
        #get pixel coordinates:
        xp,yp=annim.topixel(ra,dec)
        annim.Marker(x=xp,y=yp,
                     marker='o',mode='pixel',markersize=8, color=cname,
                     label = labname)

    #make figure
    annim.plot()
    #add legend
    plt.legend()
    #save fig
    plt.savefig(figname)
    plt.close()
    
    
#annim.plot()
# Set title for Matplotlib

def sky_plot(ra_array_lists,dec_array_lists,color_list,label_list,figname):
    """
    Make sky plots, based on code from Tom
    Inputs:
    ra_array_lists: list of RA arrays
    dec_array_lists: list of Dec array
    color_list: List of colors for plotting
    label_list: List of labels
    figname: output name for figure
    """
    #start the figure
    fig, ax = plt.subplots(figsize=(10,10))

    #add important lines to figure
    t = np.arange(0,2.0*np.pi,0.01)
    xx = np.cos(t)
    yy = np.sin(t)
    ax.plot(xx*10.0,yy*10.0,color='grey',linestyle='dotted')
    ax.plot(xx*30.0,yy*30.0,color='grey',linestyle='dotted')
    ax.plot(xx*60.0,yy*60.0,color='grey',linestyle='dotted')
    ax.plot(xx*45.0,yy*45.0,color='grey',linestyle='dotted')

    ax.plot([0,0],[0,-80],color='grey',linestyle='dotted')
    ax.plot([0,0],[0,80],color='grey',linestyle='dotted')
    ax.plot([0,80],[0,0],color='grey',linestyle='dotted')
    ax.plot([0,-90],[0,0],color='grey',linestyle='dotted')
    
    ax.plot([0,-45],[0,-45],color='grey',linestyle='dotted')
    ax.plot([0,45],[0,-45],color='grey',linestyle='dotted')
    ax.plot([0,45],[0,45],color='grey',linestyle='dotted')
    ax.plot([0,-45],[0,45],color='grey',linestyle='dotted')

    ax.set_xlim(-70,70)
    ax.set_ylim(-70,70)
    ax.axis('off')

    ax.text(11,0,r'$80^\circ$')
    ax.text(31,0,r'$60^\circ$')
    ax.text(45,0,r'$45^\circ$')
    ax.text(55,0,r'$30^\circ$')
    ax.text(1,-55,r'$12^{\rm h}$')
    ax.text(1,55,r'$0^{\rm h}$')

    ax.text(42,-38,r'$9^{\rm h}$')
    ax.text(-42,-38,r'$15^{\rm h}$')
    ax.text(42,38,r'$3^{\rm h}$')
    ax.text(-42,38,r'$21^{\rm h}$')

    #iterate through arrays and add to plot
    for ra,dec,cname,labname in zip(ra_array_lists,
                                    dec_array_lists,color_list,label_list):
        x = np.sin(ra/360.0*np.pi*2.0)*(90.0-dec)
        y = np.cos(ra/360.0*np.pi*2.0)*(90.0-dec)
        ax.scatter(x,y,marker='o',c=cname,label=labname)

    #add legend
    plt.legend()
    #save fig
    plt.savefig(figname)

def get_ra_dec(taskids,beams,master_tids,master_ra,master_dec,cbpos):
    """
    For a list of taskids and beams, get ra/dec for each combination
    Inputs:
    taskids (array): List/array of taskids, same length as beams
    beams (array): list/array of beams, same length as taskids
    master_tids (array): list/array of only taskids
    master_ra (array): list/array of RA matched to master_tids
    master_dec (array): list/array of Dec matched to master_tids
    cbpos (array): Array of CB offsets
    """
    ra_array = np.empty(len(taskids))
    dec_array = np.empty(len(taskids))
    #iterate through each taskid/beam pair
    for n, (tid, bm) in enumerate(zip(taskids,beams)):
        #get cbpos for that beam
        cb_ra_off = cbpos['ra'][bm]
        cb_dec_off = cbpos['dec'][bm]
        master_ind = np.where(master_tids == tid)[0]
        coord_field = SkyCoord(master_ra[master_ind],master_dec[master_ind],
                               frame='icrs',unit='deg')
        #have offsets, need to add properly
        #then like to check
        dec_beam = coord_field.dec.deg + cb_dec_off
        #print(np.cos(coord_field.dec.to(u.radian)))
        ra_beam = coord_field.ra.deg + cb_ra_off / np.cos(coord_field.dec.to(u.radian).value)
        coord_beam = SkyCoord(ra_beam,dec_beam,
                              frame = 'icrs', unit = 'deg')
        dra,ddec = coord_field.spherical_offsets_to(coord_beam)
        #print(cb_ra_off, cb_ra_off / np.cos(coord_field.dec.to(u.radian).value), dra.deg )
        #check how well i'm doing
        #a beam is 30', 0.5 deg. Want to be w/in 10%, so 0.05 deg
        radiff = dra.deg - cb_ra_off
        decdiff = ddec.deg - cb_dec_off
        if radiff > 0.05:
            print(('Taskid {0}, beam {1} '
                   'fails RA coordinate tolerance').format(tid,bm))
        if decdiff > 0.05:
            print(('Taskid {0}, beam {1} '
                   'fails Dec coordinate tolerance').format(tid,bm))
        ra_array[n] = coord_beam.ra.deg
        dec_array[n] = coord_beam.dec.deg

    return ra_array, dec_array

        

        
