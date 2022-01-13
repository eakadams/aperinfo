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
from functions.plots import plot_hist

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
        self.continfo = ascii.read(contfile)
        self.polinfo = ascii.read(polfile)
        self.hiinfo = ascii.read(hifile)

        """
        #remove anything that is from before survey start
        #this is code specifically for handling survey observations
        ind_survey = np.where(self.continfo['taskid'] > 190702001)[0]
        self.continfo = self.continfo[ind_survey]

        ind_survey = np.where(self.polinfo['taskid'] > 190702001)[0]
        self.polinfo = self.polinfo[ind_survey]
        """
        
        #setup beam info based on obsinfo
        beam_table_length = len(self.obsinfo) * 40
        #col_info = [ ('ObsID', 'i4'), ('beam', 'i4'), ('BeamID', 'U12'),
        #             ('cont_pass', "?"), ('cont_s_in', 'f8'),
        #             ('cont_s_out', 'f8'), ('cont_rat','f8'),
        #             ('cont_N2','f8'), ('cont_EX-2','f8')]
        #self.beaminfo = Table(data = np.full(beam_table_length, np.nan,
        #                                     dtype = col_info) )
        #boolean column defaults to True. Change that to False (have to update to Tru)
        #self.beaminfo['cont_pass'] = np.full(beam_table_length, False)
        self.beaminfo = Table()
        self.beaminfo['ObsID'] = np.empty(beam_table_length, dtype = int)
        self.beaminfo['beam'] = np.empty(beam_table_length, dtype = int)
        for n,tid in enumerate(self.obsinfo['taskID']):
            ind = n*40
            self.beaminfo['ObsID'][ind:ind+40] = np.full(40,tid)
            self.beaminfo['beam'][ind:ind+40] = np.arange(40)

        #create a "beamid" for all tables that is ObsID_beam
        #this is unique identifier for matching
        self.continfo['BeamID'] = [ f"{tid}_{b:02d}" for tid,b in
                                    self.continfo['taskid','beam'] ]
        self.beaminfo['BeamID'] = [ f"{tid}_{b:02d}" for tid, b in
                                    self.beaminfo['ObsID','beam'] ]
        self.polinfo['BeamID'] = [ f"{tid}_{b:02d}" for tid,b in
                                    self.polinfo['taskid','beam'] ]
        self.hiinfo['BeamID'] = [ f"{tid}_{b:02d}" for tid,b in
                                    self.hiinfo['Obsid','Beam'] ]
        
        #join continuum table to full beaminfo table
        #first drop duplicate columns I don't want
        self.continfo.remove_columns(['taskid','beam'])
        self.beaminfo = join(self.beaminfo, self.continfo,
                             keys = 'BeamID', join_type = 'left',
                             table_names = ['Beams', 'cont'] )

        #now add pol table
        #this is where specifying cont/pol will be important
        #since metrics are the same
        self.polinfo.remove_columns(['taskid','beam'])
        self.beaminfo = join(self.beaminfo, self.polinfo,
                             keys = 'BeamID', join_type = 'left',
                             table_names = ['cont', 'pol'])

        #and finally line table
        self.hiinfo.remove_columns(['Obsid','Beam'])
        self.beaminfo = join(self.beaminfo, self.hiinfo,
                             keys = 'BeamID', join_type = 'left')

        #note that names of columns might not alays be clearest
        #but that's not the issue of this code
        #although I may rename at some point. But it works!
   



class DR1(Beams):
    """
    Child class focused on DR1

    Attributes
    ----------

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
        Beams.__init__(self,obsfile, happilifile, contfile, hifile, polfile)

        #Limit to DR1 range
        #Do this manually by taskid

        #First exclude observations that are too late
        ind_dr1 = np.where(self.beaminfo['ObsID'] <= 200701000)[0]
        self.beaminfo = self.beaminfo[ind_dr1]

        #Then observations that are too early
        #These were not processed with 150 MHz version of pipeline
        #So all processed data is considered failed
        """
        July 2019 - 300 MHz version of pipeline
        190806345 - no processed data produced; probably due to missing set of cals at start
        190731125 - never processed
        190801041, 190801042 - never processed
        200429042 - never processed due to issues w/ ants changing for cals
        200430053-200505057 : RTC & RTD left out
        """
        ind_good_proc = np.where(
            (self.beaminfo['ObsID'] >=190807001) &
            (self.beaminfo['ObsID'] != 200429042) &
            (self.beaminfo['ObsID'] != 200430053) &
            (self.beaminfo['ObsID'] != 200501001) &
            (self.beaminfo['ObsID'] != 200501042) &
            (self.beaminfo['ObsID'] != 200502054) &
            (self.beaminfo['ObsID'] != 200503001) &
            (self.beaminfo['ObsID'] != 200503042) &
            (self.beaminfo['ObsID'] != 200505016) &
            (self.beaminfo['ObsID'] != 200505057) )[0]

        #print(len(self.beaminfo),len(ind_good_proc))
        self.beaminfo = self.beaminfo[ind_good_proc]

        #Now get just the released beams
        #Have this as a separate table because maybe I want to know
        #information about all considered
        ind_released = np.where(self.beaminfo['pass'] == 'True')[0]
        print(len(self.beaminfo),len(ind_released))

        self.released = self.beaminfo[ind_released]
        
        #Note that also need to account for skipping taskids at start w/ poor processing
        #This means should add processing information (possibly), or set manually here
        #Also want to limit to data that passed validation and hence is released.
        #So that I can focus on paper plots, will write specific plotting methods for this child object
        #But will try to call generalized functions (in plots.py) for making histograms

    def plot_cont_noise(self):
        """ 
        Plot histograms of continuum noise
        """
        fig, ax = plot_hist( self.released['s_in_cont'],
                             self.released['s_out_cont'],
                             colors = ['gray', 'black'],
                             labels = ['Inner noise', 'Outer noise'] )
        ax.legend()
        ax.set_xlabel("Noise [mJy beam$^{-1}$]")
        ax.set_ylabel("Count")
        pathname = os.path.join(figdir, 'dr1_cont_noise.pdf')
        plt.savefig(pathname)
        plt.close('all')
