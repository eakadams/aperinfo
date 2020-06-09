#Overview of Observations and Processing

from __future__ import print_function

__author__="E.A.K. Adams"

"""
Code for creating the overview
of observations and processing
Object for observations - taskid based
Object for processing - beam based

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
import info.sky_plots as sp

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
    def get_dr1_obs(self,lastind=149):
        """
        Get information for DR1, through endid
        """
        #obsinfo is already sorted
        #but need to combine with happili info
        #first, find ending index for dr1
        #didn't work, did manually, confirm it:
        print('Last taskid is {}'.format(self.obsinfo['taskID'][lastind]))
        self.dr1_obs = self.obsinfo[0:(lastind+1)]
        #check for bad data
        goodind = np.where(self.dr1_obs['quality'] == 'good')[0]
        #print(len(self.dr1_obs),len(goodind))
        #limit to good data (archived, not deleted)
        self.dr1_obs = self.dr1_obs[goodind]
        #add happili information
        (taskids, ind_dr1_obs,
         ind_happili) = np.intersect1d(self.dr1_obs['taskID'],self.happili['taskid'],
                                       return_indices=True)
        #note that there are two obs in obsinfo that hsouldnt be -- argo and early sci
        #need to check that codebut will ignore for now
        #and assume all taskids have an entry on happili
        #might be nothing due to failed processing but at least a directory exists
        #first, only keep indices that have happili entries
        self.dr1_obs = self.dr1_obs[ind_dr1_obs]
        #then add columns from happili info
        self.dr1_obs['fluxcal'] = self.happili['fluxcal'][ind_happili]
        self.dr1_obs['flux_first'] = self.happili['fluxcal_firsttaskid'][ind_happili]
        self.dr1_obs['flux_last'] = self.happili['fluxcal_lasttaskid'][ind_happili]
        self.dr1_obs['polcal'] = self.happili['polcal'][ind_happili]
        self.dr1_obs['pol_first'] = self.happili['polcal_firsttaskid'][ind_happili]
        self.dr1_obs['pol_last'] = self.happili['polcal_lasttaskid'][ind_happili]
        self.dr1_obs['apercal_version'] = self.happili['apercal_version'][ind_happili]

        #update for apercal name; need mapping
        self.dr1_obs['apercal_name'] = np.full(len(self.dr1_obs),'12345678901234567890123456')
        for i,version in enumerate(self.dr1_obs['apercal_version']):
            if version != "None":
                name = get_apercal_name(version, process = True)
                self.dr1_obs['apercal_name'][i] = name
            else:
                self.dr1_obs['apercal_name'][i] = "None"

        #check for 300 MHz processing and remove from consideration for release
        ind_good_proc = np.where(self.dr1_obs['apercal_name'] != "300 MHz")[0]
        print(len(self.dr1_obs),len(ind_good_proc))
        self.dr1_obs = self.dr1_obs[ind_good_proc]

    #make MR table for online publication
    def make_dr1_obs_mr(self):
        ascii.write(self.dr1_obs,
                    os.path.join(tabledir,'dr1_obs_mr.txt'),
                    format='cds')

    #make LATEX test table for paper
    def make_dr1_obs_paper_table(self):
        ascii.write(self.dr1_obs['taskID','name','field_ra','field_dec',
                             'fluxcal','flux_first','flux_last',
                             'polcal','pol_first','pol_last',
                             'apercal_name'][0:20],
                    os.path.join(tabledir,'dr1_obs_table_paper.txt'),
                    format='latex',
                    overwrite=True)

    #csv table for team use
    def make_dr1_obs_csv(self):
        ascii.write(self.dr1_obs,
                    os.path.join(tabledir,'dr1_obs.csv'),
                    format='csv',
                    overwrite=True)

    def plot_apercal_dr1_obs(self):
        """
        Sky plot of Apercal processing
        """
        #have to get separate lists for each Apercal processing
        #get unique names
        apercal_names = np.unique(self.dr1_obs['apercal_name'])
        #get colors
        prop_cycle = plt.rcParams['axes.prop_cycle']
        mpcolors = prop_cycle.by_key()['color']
        colorlist = mpcolors[0:len(apercal_names)]
        #get ra and dec list
        ralist = []
        declist = []
        for name in apercal_names:
            ind = np.where(self.dr1_obs['apercal_name'] == name)[0]
            ra = self.dr1_obs['field_ra'][ind]
            dec = self.dr1_obs['field_dec'][ind]
            ralist.append(ra)
            declist.append(dec)
        #make the plot
        sp.sky_plot_kapteyn(ralist,
                            declist,
                            colorlist,
                            apercal_names,
                            os.path.join(figdir,'apercal_processing_dr1_obs.pdf'))

#define class that is beam-based for processed data
#inherits from ObsCat so that I can keep that work
class ProcCat(ObsCat):
    def __init__(self,
                 obsfile = os.path.join(filedir,'obsatdb.csv'),
                 happilifile = os.path.join(filedir,'happili.csv'),
                 validfile = os.path.join(filedir,'combined_valid.csv')):
        """
        Initalize ProcCat object
        This has information about processing (per beam)
        Is child class from ObsCat so inherits all that information
        Inputs:
        - obsfile (str): File with observation information (taskid, ra/dec). Default to standard w/in package
        """
        #add additional initializion
        ObsCat.__init__(self,obsfile,happilifile)
        self.valid = ascii.read(validfile)

    #can I define a new get dr1_obs that does what i want here?
    #get DR1_OBS data
    def get_dr1_proc(self):
        """
        Get information for DR1 processing
        """
        #first get obs info
        self.get_dr1_obs()
        #now get the corresponding beam information
        #only want infomration for taskids that are in dr1_obs
        #and which pass validation
        #validation is only on continuum currently
        #need to carry polarization and line flags
        #add a field name to valid; want this info
        self.valid['Field'] = np.full(len(self.valid),'X0000+9999')
        #initialize empty array for holding indices
        array_passind = np.empty(0,dtype=int)
        #iterate through obs task ids
        for i,taskid in enumerate(self.dr1_obs['taskID']):
            passind = np.where((self.valid['taskid'] == taskid ) &
                               (self.valid['cont_pass'] == 'True') )[0]
            if len(passind) > 0:
                #if there are things that pass on to release, append them
                array_passind = np.append(array_passind,passind)

            #also update field name
            taskind = np.where(self.valid['taskid'] == taskid)[0]
            if len(taskind) > 0:
                self.valid['Field'][taskind] = self.dr1_obs['name'][i]

        #now I want to keep just beams that pass validation and are in release
        self.dr1_proc = self.valid[array_passind]

        print(len(self.valid),len(self.dr1_proc))

    #make LATEX test table for paper
    def make_dr1_proc_paper_table(self):
        ascii.write(self.dr1_proc['taskid','beam','Field',
                                  'pol_V_pass','pol_QU_pass',
                                  'HI_all_good','HI_c2_good'][0:20],
                    os.path.join(tabledir,'dr1_proc_table_paper.txt'),
                    format='latex',
                    overwrite=True)

    #csv table for team use
    def make_dr1_proc_csv(self):
        ascii.write(self.dr1_proc,
                    os.path.join(tabledir,'dr1_proc.csv'),
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
