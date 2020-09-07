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
import pandas as pd

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
#at a differnt level this time, so actually in this dir
aperinfodir = this_dir
filedir = os.path.join(aperinfodir,"files")
tabledir = os.path.join(aperinfodir,"tables")
figdir = os.path.join(aperinfodir,"figures")
#print(filedir)

#get mpl colors
prop_cycle = plt.rcParams['axes.prop_cycle']
mpcolors = prop_cycle.by_key()['color']


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
        #merge happili info into obsinfo
        #get indices
        (taskids, ind_obs,
         ind_happili) = np.intersect1d(self.obsinfo['taskID'],self.happili['taskid'],
                                       return_indices=True)
        #create columns for happili info
        self.obsinfo['fluxcal'] = np.full(len(self.obsinfo),'3C???')
        self.obsinfo['flux_first'] = np.full(len(self.obsinfo),'YYMMDDNNN')
        self.obsinfo['flux_last'] = np.full(len(self.obsinfo),'YYMMDDNNN')
        self.obsinfo['polcal'] = np.full(len(self.obsinfo),'3C???')
        self.obsinfo['pol_first'] = np.full(len(self.obsinfo),'YYMMDDNNN')
        self.obsinfo['pol_last'] = np.full(len(self.obsinfo),'YYMMDDNNN')
        self.obsinfo['apercal_version'] = np.full(len(self.obsinfo),'v0.0-???-githash-')
        #then add columns from happili info
        #first have to create colum
        self.obsinfo['fluxcal'][ind_obs] = self.happili['fluxcal'][ind_happili]
        self.obsinfo['flux_first'][ind_obs] = self.happili['fluxcal_firsttaskid'][ind_happili]
        self.obsinfo['flux_last'][ind_obs] = self.happili['fluxcal_lasttaskid'][ind_happili]
        self.obsinfo['polcal'][ind_obs] = self.happili['polcal'][ind_happili]
        self.obsinfo['pol_first'][ind_obs] = self.happili['polcal_firsttaskid'][ind_happili]
        self.obsinfo['pol_last'][ind_obs] = self.happili['polcal_lasttaskid'][ind_happili]
        self.obsinfo['apercal_version'][ind_obs] = self.happili['apercal_version'][ind_happili]

        #update for apercal name; need mapping
        self.obsinfo['apercal_name'] = np.empty(len(self.obsinfo),dtype='U30')
        for i,version in enumerate(self.obsinfo['apercal_version']):
            #print(version)
            if version == "v0.0-???-githash-":
                self.obsinfo['apercal_name'][i] = "Not processed"
            elif version != "None":
                name = get_apercal_name(version, process = True)
                self.obsinfo['apercal_name'][i] = name
            else:
                self.obsinfo['apercal_name'][i] = "None"

        #add tex name for later
        self.obsinfo['apercal_tex_name'] = self.obsinfo['apercal_name']
        for i,name in enumerate(self.obsinfo['apercal_tex_name']):
            self.obsinfo['apercal_tex_name'][i] = name.replace("_","\_")
            

        #add observational notes
        self.obs_notes = ascii.read(os.path.join(filedir,'obs_notes.csv'))
        self.obsinfo['reobserve'] = np.full(len(self.obsinfo),'N')
        self.obsinfo['reprocess'] = np.full(len(self.obsinfo),'N')
        #object for arbitrary string length
        self.obsinfo['note'] = np.empty(len(self.obsinfo),dtype='object')
        (taskids, ind_obs,
         ind_note) = np.intersect1d(self.obsinfo['taskID'],self.obs_notes['taskid'],
                                       return_indices=True)
        self.obsinfo['reobserve'][ind_obs] = self.obs_notes['reobserve'][ind_note]
        self.obsinfo['reprocess'][ind_obs] = self.obs_notes['reprocess'][ind_note]
        self.obsinfo['note'][ind_obs] = self.obs_notes['note'][ind_note]


    #get obs for data release
    def get_data_release_obs(self,firstind=0,lastind=148):
        """
        Get obs specified by firstind/lastind
        """
        print('First taskid is {}'.format(self.obsinfo['taskID'][firstind]))
        print('Last taskid is {}'.format(self.obsinfo['taskID'][lastind]))
        data_release = self.obsinfo[firstind:(lastind+1)]
        #check for bad data
        goodind = np.where(data_release['quality'] == 'good')[0]
        #print(len(self.dr1_obs),len(goodind))
        #limit to good data (archived, not deleted)
        data_release = data_release[goodind]
        #check for ahppili info
        (taskids, ind_dr1_obs,
         ind_happili) = np.intersect1d(data_release['taskID'],self.happili['taskid'],
                                       return_indices=True)
        #note that there are two obs in obsinfo that hsouldnt be -- argo and early sci
        #need to check that codebut will ignore for now
        #and assume all taskids have an entry on happili
        #might be nothing due to failed processing but at least a directory exists
        #first, only keep indices that have happili entries
        data_release = data_release[ind_dr1_obs]
        return data_release

    #get dr1 data
    def get_dr1_obs(self):
        self.dr1_obs = self.get_data_release_obs(firstind=0,lastind=148)

    #get dr1 plus
    def get_dr1plus_obs(self):
        self.dr1_plus_obs = self.get_data_release_obs(firstind=0,lastind=180)

       
#    #get DR1 data
#    def get_dr1_obs(self,lastind=148):
#        """
#        Get information for DR1, through endid
#        """
#        #obsinfo is already sorted
#        #and now have added happili info there already
#        #first, find ending index for dr1
#        #didn't work, did manually, confirm it:
#        print('Last taskid is {}'.format(self.obsinfo['taskID'][lastind]))
#        self.dr1_obs = self.obsinfo[0:(lastind+1)]
#        #check for bad data
#        goodind = np.where(self.dr1_obs['quality'] == 'good')[0]
#        #print(len(self.dr1_obs),len(goodind))
#        #limit to good data (archived, not deleted)
#        self.dr1_obs = self.dr1_obs[goodind]
#        #check for ahppili info
#        (taskids, ind_dr1_obs,
#         ind_happili) = np.intersect1d(self.dr1_obs['taskID'],self.happili['taskid'],
#                                       return_indices=True)
#        #note that there are two obs in obsinfo that hsouldnt be -- argo and early sci
#        #need to check that codebut will ignore for now
#        #and assume all taskids have an entry on happili
#        #might be nothing due to failed processing but at least a directory exists
#        #first, only keep indices that have happili entries
#        self.dr1_obs = self.dr1_obs[ind_dr1_obs]
#
#        """
#        This is only for processed data, still releasing for raw 
#        observational data
#        #check for 300 MHz processing and remove from consideration for release
#        ind_good_proc = np.where(self.dr1_obs['apercal_name'] != "300 MHz")[0]
#        print(len(self.dr1_obs),len(ind_good_proc))
#        self.dr1_obs = self.dr1_obs[ind_good_proc]
#        """


        
    #make MR table for online publication
    def make_dr1_obs_mr(self):
        ascii.write(self.dr1_obs,
                    os.path.join(tabledir,'dr1_obs_mr.txt'),
                    format='cds',overwrite=True)

    #make LATEX test table for paper
    def make_dr1_obs_paper_table(self):
        #set column names to be tex & user friendly
        col_names = ['ObsID','Name','RA','Dec','Fluxcal','flux\_first','flux\_last',
                     'Polcal','pol\_first','pol\_last','Apercal\_name','Apercal\_version']
        
        ascii.write(self.dr1_obs['taskID','name','field_ra','field_dec',
                             'fluxcal','flux_first','flux_last',
                             'polcal','pol_first','pol_last',
                                 'apercal_tex_name','apercal_version'][0:30],
                    os.path.join(tabledir,'dr1_obs_table_paper.txt'),
                    format='latex',
                    overwrite=True,
                    names=col_names,
                    col_align=len(col_names)*'l',
                    latexdict = {'header_start': "\hline \hline",
                                 'header_end': "\hline",
                                 'data_end': "\hline",
                                 'caption': "Summary of released survey observation",
                                 'preamble': ["\centering","\label{tab:obs}"]}
        )

    #make obs_notes table
    def make_dr1_obs_notes_table(self):
        """
        Make a latex formatted table of obsnotes for paper.
        Example 30 rows; part of obs table in machine-readable version
        """
        ascii.write(self.dr1_obs['taskID','name','note'][0:30],
                    os.path.join(tabledir,'dr1_obs_notes.txt'),
                    format='latex',
                    overwrite=True,
                    names=['ObsID','Name','Notes'],
                    col_align='llp{12cm}',
                    latexdict = {'header_start': "\hline \hline",
                                 'header_end': "\hline",
                                 'data_end': "\hline",
                                 'caption': "Example of observation notes",
                                 'preamble': ["\centering","\label{tab:obsnotes}"],
                                 'tabletype': "table*"}
                    )
        

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
        #add a name for medium-deep
        apercal_names = np.append(apercal_names,'AMES')
        #get colors
        prop_cycle = plt.rcParams['axes.prop_cycle']
        mpcolors = prop_cycle.by_key()['color']
        colorlist = mpcolors[0:len(apercal_names)]
        #get ra and dec list; need to first separate medium-deep pointings
        ind_ames = [i for i, s in enumerate(self.dr1_obs['name']) if 'M' in s]
        ind_awes = [i for i, s in enumerate(self.dr1_obs['name']) if 'S' in s]
        self.dr1_ames =  self.dr1_obs[ind_ames]
        self.dr1_awes =  self.dr1_obs[ind_awes]


        #also want to find only those that have duplicates, e.g., multiple observations
        s = pd.Series(self.dr1_ames['name'])
        dup = s[s.duplicated()]
        repeated_fields = np.unique(dup)


        #need to find part dr1_ames that contains repeated_fields above
        #match two string arrays...
        #almost could get info from pd dup object but it doesn't count first occurence
        ind_repeats = []
        for field in repeated_fields:
            ind = [i for i, s in enumerate(self.dr1_ames['name']) if field in s]
            ind_repeats.append(ind)

        repeat_array = np.sort(np.hstack(ind_repeats))
        #np.sort(ind_repeats)
        print(repeat_array)

        self.dr1_repeated_ames = self.dr1_ames[repeat_array]
        
        ralist = []
        declist = []
        for name in apercal_names[0:-1]:
            ind = np.where(self.dr1_obs['apercal_name'] == name)[0]
            ra = self.dr1_obs['field_ra'][ind]
            dec = self.dr1_obs['field_dec'][ind]
            ralist.append(ra)
            declist.append(dec)

        #add repeated ames to end
        ra = self.dr1_repeated_ames['field_ra']
        dec = self.dr1_repeated_ames['field_dec']
        ralist.append(ra)
        declist.append(dec)

        
        #make the plots
        #want to separate medium-deep points so can plot separately
        #all sky plot
        sp.sky_plot_kapteyn(ralist,
                            declist,
                            colorlist,
                            apercal_names,
                            os.path.join(figdir,'apercal_processing_dr1_obs.pdf'),
                            survey_pointings = os.path.join(filedir,'all_pointings.v7.18jun20.txt'))

        #need to add a separate medium-deep plot
        #first sort by field name
        field_name = np.flip(np.sort(np.unique(self.dr1_repeated_ames['name'])))

        print(len(self.dr1_repeated_ames))

        #then create plot coordinates for everything
        #base on field name, want to be in same row
        plot_x = np.full(len(self.dr1_repeated_ames),-10)
        plot_y = np.full(len(self.dr1_repeated_ames),-10)

        for i, field in enumerate(field_name):
            ames_ind = np.where(self.dr1_repeated_ames['name'] == field)[0]
            #all fields in same row, same y coord
            #offset by 1 for ease of plotting
            plot_y[ames_ind] = i+1
            #for x coordinate, have to iterate through all fields
            #offset by 1 again for plot legibility
            for k in range(len(ames_ind)):
                plot_x[ames_ind[k]] = k+1

        #setup figure
        #fig, ax = plt.subplots(figsize=[6,4])
        fig = plt.figure(figsize=[6,4])
        ax = fig.add_axes([0.2, 0.15, .75, .75 ])


        #have coordinates for al fields, now have to iterate over Apercal name for colors
        #skip last one, placeholder for AMES
        for color,name in zip(colorlist[0:-1],apercal_names[0:-1]):
            plotind = np.where(self.dr1_repeated_ames['apercal_name'] == name)[0]
            ax.scatter(plot_x[plotind],plot_y[plotind],c=color,label=name)

        #plt.legend()

        ax.set_yticks(list(range(1,len(field_name)+1)))
        ax.set_yticklabels(list(field_name))

        ax.set_title('Medium-deep fields')
        ax.set_xlabel('Number of observations')

        plotname = os.path.join(figdir,'apercal_processing_dr1_ames.pdf')
        plt.savefig(plotname)
        plt.close()

        
    def plot_all_obs(self):
        """
        Make a sky plot of all observations
        Won't worry about MDS at this time, 
        just want to see sky coverage
        Will still color by Apercal processing
        Useful to see what is good processing
        and what hasn't been processed
        """
        #have to get separate lists for each Apercal processing
        #get unique names
        apercal_names = np.unique(self.obsinfo['apercal_name'])
        #get colors
        prop_cycle = plt.rcParams['axes.prop_cycle']
        mpcolors = prop_cycle.by_key()['color']
        colorlist = mpcolors[0:len(apercal_names)]

        #get ra/dec lists for plotting
        ralist = []
        declist = []
        for name in apercal_names:
            ind = np.where(self.obsinfo['apercal_name'] == name)[0]
            ra = self.obsinfo['field_ra'][ind]
            dec = self.obsinfo['field_dec'][ind]
            ralist.append(ra)
            declist.append(dec)

        #make the plots
        #want to separate medium-deep points so can plot separately
        #all sky plot
        sp.sky_plot_kapteyn(ralist,
                            declist,
                            colorlist,
                            apercal_names,
                            os.path.join(figdir,'apercal_all_obs.pdf'),
                            survey_pointings = os.path.join(filedir,'all_pointings.v7.18jun20.txt'))
        
        plt.close()
        

    def plot_obs_by_month(self):
        """
        Sky plot of observations, color coded by time
        Take everything w/in DR1 as one color, then plot by month after that
        Don't use a full alpha so can see mds a bit - still need to do this
        Manually set months and update as needed:
        Feb, Mar, Apr, May, June
        Know everything has been processed already
        """
        month_list = ['DR1','02','03','04','05','06']
        colorlist = mpcolors[0:len(month_list)]

        ralist = []
        declist = []
        for month in month_list:
            if month == 'DR1':
                ra = self.dr1_obs['field_ra']
                dec = self.dr1_obs['field_dec']
            else:
                yymm = '20'+month
                str_taskid = self.obsinfo['taskID'].astype(str)
                ind = np.flatnonzero(np.core.defchararray.find(
                    str_taskid,
                    yymm,
                    end=4)!=-1)
                ra = self.obsinfo['field_ra'][ind]
                dec = self.obsinfo['field_dec'][ind]
            ralist.append(ra)
            declist.append(dec)

        #all sky plot
        sp.sky_plot_kapteyn(ralist,
                            declist,
                            colorlist,
                            month_list,
                            os.path.join(figdir,'apercal_obs_by_month.pdf'),
                            survey_pointings = os.path.join(filedir,'all_pointings.v7.18jun20.txt'))
        
        plt.close()

        
    def plot_dr1_context(self):
        """
        Make a plot of sky  coverage that is color-coded by DR1 and then by  month
        """
        ra_dr1 = self.dr1_obs['ra']
        dec_dr1 = self.dr1_obs['dec']

        #get feb, mar, apr, may
        #june is incomplete processing
        #plot w/ dr1 last to deal with feb issue (only 1-3 fields)
        

#define class that is beam-based for processed data
#inherits from ObsCat so that I can keep that work
class ProcCat(ObsCat):
    def __init__(self,
                 obsfile = os.path.join(filedir,'obsatdb.csv'),
                 happilifile = os.path.join(filedir,'happili.csv'),
                 validfile = os.path.join(filedir,'combined_valid.csv'),
                 beampos = os.path.join(filedir,'cb_offsets.txt')):
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
        #and beam position
        self.cb_pos = ascii.read(beampos)

    #can I define a new get dr1_obs that does what i want here?
    #get DR1_OBS data
    def get_dr1_proc(self):
        """
        Get information for DR1 processing
        """
        #first get obs info
        self.get_dr1_obs()
        #check for 300 MHz processing and remove from consideration for release
        #also for data sets that weren't processed
        #and any other where I have notes:
        """
        190806345 - no processed data produced; probably due to missing set of cals at start
        191024043 - RT2 delays bad; need to flag and reprocess
        191024044 - RT2 delays bad; need to flag and reprocess
        191209025 - RT9 off source in middle of observation; flag and reprocess
        191227014 - RT9 not on source; flag and reprocess
        191230041 - RT9 not on source for ~first hour; 
                    RT8 had large self-cal amp solutions after that. Flag & reprocess
        191031242 - bright stripes --> bad time at start to be flagged. Also RT3
        191209026 - bright sources on edge not in mask. Need updated pipeline version
        """
        ind_good_proc = np.where((self.dr1_obs['apercal_name'] != "Apercal_300") &
                                 (self.dr1_obs['apercal_name'] != "None") &
                                 (self.dr1_obs['taskID'] != 190806345) &
                                 (self.dr1_obs['taskID'] != 191024043) &
                                 (self.dr1_obs['taskID'] != 191024044) &
                                 (self.dr1_obs['taskID'] != 191209025) &
                                 (self.dr1_obs['taskID'] != 191227014) &
                                 (self.dr1_obs['taskID'] != 191230041) &
                                 (self.dr1_obs['taskID'] != 191031242) &
                                 (self.dr1_obs['taskID'] != 191209026))[0]
        print(len(self.dr1_obs),len(ind_good_proc))
        self.dr1_obs = self.dr1_obs[ind_good_proc]
        #now get the corresponding beam information
        #only want infomration for taskids that are in dr1_obs
        #and which pass validation
        #validation is only on continuum currently
        #need to carry polarization and line flags
        #add a field name to valid; want this info
        self.valid['Field'] = np.full(len(self.valid),'X0000+9999')
        #add a N_pass to obs table
        self.dr1_obs['N_pass'] = np.full(len(self.dr1_obs),0)
        #initialize empty array for holding indices
        array_passind = np.empty(0,dtype=int)
        #iterate through obs task ids
        for i,taskid in enumerate(self.dr1_obs['taskID']):
            passind = np.where((self.valid['taskid'] == taskid ) &
                               (self.valid['cont_pass'] == 'True') )[0]
            #add number that pass
            self.dr1_obs['N_pass'][i] = len(passind)
            if len(passind) > 0:
                #if there are things that pass on to release, append them
                array_passind = np.append(array_passind,passind)

            #also update field name
            taskind = np.where(self.valid['taskid'] == taskid)[0]
            if len(taskind) > 0:
                self.valid['Field'][taskind] = self.dr1_obs['name'][i]

        #now I want to keep just beams that pass validation and are in release
        self.dr1_proc = self.valid[array_passind]

        #add ra and dec
        ra_array, dec_array = sp.get_ra_dec(self.dr1_proc['taskid'],
                                            self.dr1_proc['beam'],
                                            self.dr1_obs['taskID'],
                                            self.dr1_obs['field_ra'],
                                            self.dr1_obs['field_dec'],
                                            self.cb_pos)
        self.dr1_proc['ra'] = ra_array
        self.dr1_proc['dec'] = dec_array

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

    #make continuum table for paper
    def make_dr1_cont_table(self):
        #set column names to be tex & user friendly
        col_names = ['ObsID','Name','Beam','$\sigma_{in}$',
                     '$\sigma_{out}$', '$R$', 'N2', 'Ex-2']
        
        ascii.write(self.dr1_proc['taskid','Field','beam','s_in','s_out',
                                 'rat','N2','Ex-2'][0:30],
                    os.path.join(tabledir,'dr1_cont.txt'),
                    format='latex',
                    overwrite=True,
                    names=col_names,
                    col_align=len(col_names)*'l',
                    latexdict = {'header_start': "\hline \hline",
                                 'header_end': "\hline",
                                 'data_end': "\hline",
                                 'caption': "Continuum validation metrics for released data",
                                 'preamble': ["\centering","\label{tab:cont}","\small"]}
                    )

    #make pol table for paper
    def make_dr1_pol_table(self):
        #set column names to be tex & user friendly
        col_names = ['ObsID','Name','Beam','Status','V status','QU status',
                     '$\sigma_{in}$',
                     '$\sigma_{out}$', '$R$', 'N2', 'Ex-2',
                     '$FT_{max}$','$p_{in}$','P2','$b_{min}$',
                     '$Q_{beam}$','$U_{beam}$','$Q_{noise}$','$U_{noise}$']
        
        ascii.write(self.dr1_proc['taskid','Field','beam',
                                 'pol_pass','pol_V_pass','pol_QU_pass',
                                 'pol_s_in','pol_s_out',
                                 'pol_rat','pol_N2','pol_Ex-2','pol_ftmax',
                                 'pol_peak_in','pol_P2','pol_bmin','Q_bm_fg',
                                 'U_bm_fg', 'Q_st_fg','U_st_fg'][0:30],
                    os.path.join(tabledir,'dr1_pol.txt'),
                    format='latex',
                    overwrite=True,
                    names=col_names,
                    col_align=len(col_names)*'l',
                    latexdict = {'header_start': "\hline \hline",
                                 'header_end': "\hline",
                                 'data_end': "\hline",
                                 'caption': "Polarization validation metrics for released data",
                                 'preamble': ["\centering","\label{tab:pol}"]}
                    )
                    

    #make HI table for paper
    def make_dr1_hi_table(self):
    #set column names to be tex & user friendly
        col_names = ['ObsID','Name','Beam','All good', 'All good/okay',
                     'cube2', 'cube1', 'cube0','$\sigma_{c2}$',
                     '$\sigma_{c1}$', '$\sigma_{c0}$',
                     '$f_{ex,c2}$', '$f_{ex,c1}$', '$f_{ex,c0}$',
                     '$p_{0.8,c2}$','$p_{0.8,c1}$','$p_{0.8,c0}$']

        #put sigma into mJy, rather than Jy units before formatting
        self.dr1_proc['rms_c1_mJy'] = 1e3*self.dr1_proc['rms_c1']
        self.dr1_proc['rms_c2_mJy'] = 1e3*self.dr1_proc['rms_c2']
        self.dr1_proc['rms_c0_mJy'] = 1e3*self.dr1_proc['rms_c0']
        
        ascii.write(self.dr1_proc['taskid','Field','beam','HI_all_good',
                                 'HI_all_good_ok','c2','c1','c0','rms_c2_mJy',
                                 'rms_c1_mJy','rms_c0_mJy','lgfrac_c2','lgfrac_c1',
                                 'lgfrac_c0','prom_c2','prom_c1','prom_c0'][0:30],
                    os.path.join(tabledir,'dr1_line.txt'),
                    format='latex',
                    overwrite=True,
                    names=col_names,
                    col_align=len(col_names)*'l',
                    latexdict = {'header_start': "\hline \hline",
                                 'header_end': "\hline",
                                 'data_end': "\hline",
                                 'caption': "Line validation metrics for released data",
                                 'preamble': ["\centering","\label{tab:hi}"]},
                    formats = {'$\sigma_{c2}$': '4.2f', '$\sigma_{c1}$': '4.2f',
                               '$\sigma_{c0}$': '4.2f', '$f_{ex,c2}$': '5.2f',
                               '$f_{ex,c1}$': '5.2f', '$f_{ex,c0}$': '5.2f',
                               '$p_{0.8,c2}$': '4.2f', '$p_{0.8,c1}$': '4.2f',
                               '$p_{0.8,c0}$': '4.2f'}
                    )

    

    def plot_dr1_proc(self):
        """
        Make sky plots of the processed data
        Will start by showing overall sky coverage, valid cont data
        Then for released data want plots that show good/okay/bad pol & line
        """
        prop_cycle = plt.rcParams['axes.prop_cycle']
        mpcolors = prop_cycle.by_key()['color']
        colorlist = [mpcolors[0]]

        plot_processed_data([self.dr1_proc['ra']],
                            [self.dr1_proc['dec']],
                            colorlist,
                            ['Released beams'],
                            self.dr1_obs['taskID','field_ra','field_dec'],
                            'dr1_proc')

        #plot Stokes V quality
        quality = ['Pass','Fail']
        colorlist = mpcolors[0:len(quality)]
        ralist = []
        declist = []
        for qual in quality:
            if qual == 'Pass':
                ind = np.where(self.dr1_proc['pol_V_pass'] == 'True')[0]
            if qual == 'Fail':
                ind = np.where(self.dr1_proc['pol_V_pass'] == 'False')[0]
            ra = self.dr1_proc['ra'][ind]
            dec = self.dr1_proc['dec'][ind]
            ralist.append(ra)
            declist.append(dec)

        plot_processed_data(ralist,declist,colorlist,quality,
                            self.dr1_obs['taskID','field_ra','field_dec'],'dr1_proc_V')

        
        #plot Stokes QU quality
        quality = ['Pass','Fail']
        colorlist = mpcolors[0:len(quality)]
        ralist = []
        declist = []
        for qual in quality:
            if qual == 'Pass':
                ind = np.where(self.dr1_proc['pol_QU_pass'] == 'True')[0]
            if qual == 'Fail':
                ind = np.where(self.dr1_proc['pol_QU_pass'] == 'False')[0]
            ra = self.dr1_proc['ra'][ind]
            dec = self.dr1_proc['dec'][ind]
            ralist.append(ra)
            declist.append(dec)

        plot_processed_data(ralist,declist,colorlist,quality,
                            self.dr1_obs['taskID','field_ra','field_dec'],'dr1_proc_QU')

        
        #plot HI data quality
        #iterate through cubes
        cubes = ['c0','c1','c2']
        quality = ['Good','Bad','Okay']
        colorlist = mpcolors[0:len(quality)]
        for cube in cubes:
            ralist = []
            declist = []
            for qual in quality:
                if qual == 'Good':
                    ind = np.where(self.dr1_proc[cube] == 'G')[0]
                if qual == 'Okay':
                    ind = np.where(self.dr1_proc[cube] == 'O')[0]
                if qual == 'Bad':
                    ind = np.where(self.dr1_proc[cube] == 'B')[0]
                ra = self.dr1_proc['ra'][ind]
                dec = self.dr1_proc['dec'][ind]
                ralist.append(ra)
                declist.append(dec)

            plot_processed_data(ralist,declist,colorlist,quality,
                                self.dr1_obs['taskID','field_ra','field_dec'],'dr1_proc_HI'+cube)


    #explore DR1
    def explore_dr1(self):
        """
        Explore the DR1 options
        Report various numbers / statistics
        """
        #first get total fraction of beams that pass
        nbeam_total = len(self.dr1_obs)*40
        nbeam_pass = len(self.dr1_proc)
        print(("{0} beams of {1} possible passed, "
               "for a total of {2:2.0f}%").format(nbeam_pass,nbeam_total,
                                             nbeam_pass/nbeam_total*100))

        #now update to remove observations with known issues
        (bad_taskids, ind_obs,
         ind_notes) = np.intersect1d(self.dr1_obs['taskID'],self.obs_notes['taskid'],
                                     return_indices=True)
        good_taskids = np.setxor1d(bad_taskids,self.dr1_obs['taskID'])
        mask = np.isin(self.dr1_proc['taskid'],good_taskids)
        #print(len(mask))
        #now find how many beams considered here
        nbeam_good = len(np.where(mask == True)[0])
        nbeam_poss = len(good_taskids)*40
        print(("When considering only good observations,"
               "{0} beams of {1} possible passed, "
               "for a total of {2:2.0f}%").format(nbeam_good,nbeam_poss,
                                             nbeam_good/nbeam_poss*100))
        #
        #c = np.setdiff1d(np.union1d(a, b), np.intersect1d(a, b))

        #add further constraints here - take data in certain time ranges, for example
        #based on system improvements / processing

        #search for best taskids
        #want taskids where fewer than 4 beams fail continuum validation
        #should add this information to dr1_obs table; easiest way to do things
        ind_best = np.where(self.dr1_obs['N_pass'] >= 32)[0]
        print(("There are {0} fields that have 32 or more beams "
               "pass continuum validation").format(len(ind_best)))
        #how do these fields do on HI / polarization validation?
        best_task = self.dr1_obs['taskID'][ind_best]
        print(("The best taskids are {0}").format(best_task))
        print(("The best fields are {0}").format(self.dr1_obs['name'][ind_best]))
        mask_best = np.isin(self.dr1_proc['taskid'],best_task)
        self.dr1_best_proc = self.dr1_proc[mask_best]
        pass_cont = np.where(self.dr1_best_proc['cont_pass'] == 'True')[0]
        print(("There are {0} possible beams "
               "in best processed fields. "
               "{1} have images. {2} "
               "pass continuum validation").format(40*len(best_task),len(self.dr1_best_proc),len(pass_cont)))
        good_HI_ind = np.where(self.dr1_best_proc['HI_c2_good_ok'] == 'True')[0]
        print(("{0} beams have good/ok HI "
               "in cube 2").format(len(good_HI_ind)))
        pass_pol = np.where(self.dr1_best_proc['pol_pass'] == 'True')[0]
        print(("{0} pass polarization validation").format(len(pass_pol)))
        
            
            
        
    


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


def plot_processed_data(ralist,declist,colorlist,labellist,obsdata,figbase):
    """
    A helper  function to produce sky plots for processed data
    Takes ralist, declist, colorlist, labellist as inputs for plotting
    also takes an array/list of info for plotting all obs for reference
    And finally a base name for figure, updated for spring & fall views
    Produces fully sky, spring and fall views
    """
    sp.sky_plot_kapteyn(ralist,
                        declist,
                        colorlist,
                        labellist,
                        os.path.join(figdir,figbase+'.pdf'),
                        mode='beam',
                        obs = obsdata
    )

    sp.sky_plot_kapteyn(ralist,
                        declist,
                        colorlist,
                        labellist,
                        os.path.join(figdir,figbase+'_spring.pdf'),
                        mode='beam',
                        obs = obsdata,
                        sky='spring'
    )

    sp.sky_plot_kapteyn(ralist,
                        declist,
                        colorlist,
                        labellist,
                        os.path.join(figdir,figbase+'_fall.pdf'),
                        mode='beam',
                        obs = obsdata,
                        sky='fall'
    )
