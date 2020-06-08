#Object for Observations

from __future__ import print_function

__author__="E.A.K. Adams"

"""
Object for observations
(ObsId, object oriented).
Used other code to be able
to make sky plots / tables

Try a class-based approach where I
create a class (centered around an
astropy Table) 
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
#at a differnt level this time, so actually in this dir
aperinfodir = this_dir
filedir = os.path.join(aperinfodir,"files")
tabledir = os.path.join(aperinfodir,"tables")
figdir = os.path.join(aperinfodir,"figures")
#print(filedir)


#define observation object
class ObsCat(object):
    def __init__(self,
                 obsfile = os.path.join(filedir,'obsatdb.csv'),
                 happilifile = os.path.join(filedir,'happili.csv')):
        """
        Initalize ObsCat object
        This has information about observaitons
        Inputs:
        - obsfile (str): File with observation information (taskid, ra/dec). Default to standard w/in package
        """
        #obsinfo as supplementary information
        self.obsinfo = ascii.read(obsfile)
        self.happili = ascii.read(happilifile)

    #get DR1 data
    def get_dr1(self,lastind=149):
        """
        Get information for DR1, through endid
        """
        #obsinfo is already sorted
        #but need to combine with happili info
        #first, find ending index for dr1
        #didn't work, did manually, confirm it:
        print('Last taskid is {}'.format(self.obsinfo['taskID'][lastind]))
        self.dr1 = self.obsinfo[0:(lastind+1)]
        #check for bad data
        goodind = np.where(self.dr1['quality'] == 'good')[0]
        #print(len(self.dr1),len(goodind))
        #limit to good data (archived, not deleted)
        self.dr1 = self.dr1[goodind]
        #add happili information
        (taskids, ind_dr1,
         ind_happili) = np.intersect1d(self.dr1['taskID'],self.happili['taskid'],
                                       return_indices=True)
        #note that there are two obs in obsinfo that hsouldnt be -- argo and early sci
        #need to check that codebut will ignore for now
        #and assume all taskids have an entry on happili
        #might be nothing due to failed processing but at least a directory exists
        #first, only keep indices that have happili entries
        self.dr1 = self.dr1[ind_dr1]
        #then add columns from happili info
        self.dr1['fluxcal'] = self.happili['fluxcal'][ind_happili]
        self.dr1['flux_first'] = self.happili['fluxcal_firsttaskid'][ind_happili]
        self.dr1['flux_last'] = self.happili['fluxcal_lasttaskid'][ind_happili]
        self.dr1['polcal'] = self.happili['polcal'][ind_happili]
        self.dr1['pol_first'] = self.happili['polcal_firsttaskid'][ind_happili]
        self.dr1['pol_last'] = self.happili['polcal_lasttaskid'][ind_happili]
        self.dr1['apercal_version'] = self.happili['apercal_version'][ind_happili]

        #update for apercal name; need mapping
        self.dr1['apercal_name'] = np.full(len(self.dr1),'12345678901234567890123456')
        for i,version in enumerate(self.dr1['apercal_version']):
            if version != "None":
                name = get_apercal_name(version, process = True)
                self.dr1['apercal_name'][i] = name
            else:
                self.dr1['apercal_name'][i] = "None"

    #make MR table for online publication
    def make_dr1_mr(self):
        ascii.write(self.dr1,
                    os.path.join(tabledir,'obs_mr.txt'),
                    format='cds')

    #make LATEX test table for paper
    def make_dr1_paper_table(self):
        ascii.write(self.dr1['taskID','name','field_ra','field_dec',
                             'fluxcal','flux_first','flux_last',
                             'polcal','pol_first','pol_last',
                             'apercal_name'][0:20],
                    os.path.join(tabledir,'obs_table_paper.txt'),
                    format='latex',
                    overwrite=True)

    #csv table for team use
    def make_dr1_obs_csv(self):
        ascii.write(self.dr1,
                    os.path.join(tabledir,'obs_dr1.csv'),
                    format='csv',
                    overwrite=True)
                    


def get_apercal_name(version,process=True):
    """
    Helper function that takes an apercal version number and 
    returns a friendly name
    Input:
        version : apercal version name
        process: optional keyword, to set names to focus on processing
    Output
        name: friendly name
    """
    #read in master file w/ information
    apercal_key = ascii.read(os.path.join(filedir,"apercal_naming.csv"))
    vind = np.where(apercal_key['Version'] == version)[0][0]
    name = apercal_key['Name'][vind]
    #redo some naming on the fly, to focus on processing
    #if proces = True
    if process is True:
        if name == "150 MHz and continuum chunks":
            name = "150 MHz"
        if name == "Add last continuum chunk image":
            name = "Bandpass phase flagging"
    return name
