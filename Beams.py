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
import tol_colors as tc
import astropy.units as u
from functions.utilities import get_beam_ra_dec
from astropy.coordinates import SkyCoord
from functions.utilities import get_med_onesig

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
#at a differnt level this time, so actually in this dir
aperinfodir = this_dir
filedir = os.path.join(aperinfodir,"files")
tabledir = os.path.join(aperinfodir,"tables")
figdir = os.path.join(aperinfodir,"figures")


#set up colors - use Paul Tol's  color blind
plt.rc('axes', prop_cycle=plt.cycler('color', list(tc.tol_cset('bright'))))
plt.cm.register_cmap('YlOrBr', tc.tol_cmap('YlOrBr'))
plt.rc('image', cmap='YlOrBr')
prop_cycle = plt.rcParams['axes.prop_cycle']
mpcolors = prop_cycle.by_key()['color']



#set fonts and latex
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

plt.rc('font', size = 20)

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
        self.beaminfo['Field'] = np.empty(beam_table_length, dtype = "S10")
        self.beaminfo['Field_RA'] = np.empty(beam_table_length)
        self.beaminfo['Field_Dec'] = np.empty(beam_table_length)
        for n,(tid,f,ra,dec) in \
            enumerate(self.obsinfo['taskID','name','field_ra','field_dec']):
            ind = n*40
            self.beaminfo['ObsID'][ind:ind+40] = np.full(40,tid)
            self.beaminfo['beam'][ind:ind+40] = np.arange(40)
            self.beaminfo['Field'][ind:ind+40] = np.full(40, f)
            self.beaminfo['Field_RA'][ind:ind+40] = np.full(40,ra)
            self.beaminfo['Field_Dec'][ind:ind+40] = np.full(40,dec)

        #add RA, Dec column to everything
        #use a helper function to get values from Field RA,Dec plus beam number
        ra, dec = get_beam_ra_dec(self.beaminfo['Field_RA'],
                                  self.beaminfo['Field_Dec'],
                                  self.beaminfo['beam'])
        self.beaminfo['RA'] = ra
        self.beaminfo['Dec'] = dec

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

    def find_radar(self, factor = 1.1):
        """ 
        Find taskids that are strongly impacted by military radar
        Use HI validation to do this, testing noise in cube2 against cube0
        Generally, cube2 should have lower noise than cube0 but 
        the radar impacts cube2 more strongly.
        use a factor to make clear things are worse.
        """
        #first fill values because masking causes me issues:
        self.beaminfo['rms_c0'].fill_value = np.nan
        self.beaminfo['rms_c2'].fill_value = np.nan

        #for reference get a list of all TIDs that have HI valid info
        ind_hi = np.where(self.beaminfo['rms_c0'].filled() > 0)[0]
        hi_valid_tids = np.unique(self.beaminfo['ObsID'][ind_hi])

        #find all beams where c2 rms greater than c0 rms
        ind_highc2 = np.where( self.beaminfo['rms_c2'].filled() > factor * self.beaminfo['rms_c0'].filled() )[0]

        print(len(ind_highc2), len(ind_hi))

        #get tids and count
        poss_radar_tids, count_bad = np.unique(
            self.beaminfo['ObsID'][ind_highc2], return_counts = True)
        
        #check for more than 10 beams showing up, then keep/report tid
        ind_radar = np.where( count_bad >= 30)[0]
        radar_tids = poss_radar_tids[ind_radar]

        return radar_tids

        
        

    def find_beam(self, ra, dec, radius = 30*u.arcmin):
        """ 
        Find beams within radius of position

        Parameters:
        -----------
        ra : float and/or Quantity
            R.A. If a float, assumed in degrees
        dec : float and/or Quantity
            Decl. If a float, assumed in degrees
        radius : float and/or Quantity
            Radius to search for beam centers. If a float, assumed in arcmins
            Default is 30'
        """
       
        #check for units
        if not isinstance(ra, u.Quantity):
            ra = ra * u.deg
        if not isinstance(dec, u.Quantity):
            dec = dec * u.deg
        if not isinstance(radius, u.Quantity):
            radius = radius * u.arcmin

        #set up coords
        all_coords = SkyCoord(self.beaminfo['RA'], self.beaminfo['Dec'],
                              frame='icrs',unit='deg')

        target_coord = SkyCoord(ra, dec, frame='icrs')

        separation = all_coords.separation(target_coord)

        idx_beams = np.where(separation < radius)[0]

        if len(idx_beams) > 0:
            beams = self.beaminfo['BeamID'][idx_beams]
            print(f"Beams within a radius of {radius} are {beams}")
        else:
            beams = None

        return beams, idx_beams

    

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
        fig, ax = plot_hist( self.released['s_out_cont'],
                             self.released['s_in_cont'],
                             colors = ['black', 'gray'],
                             labels = ['Outer noise', 'Inner noise'],
                             alpha = [1.0, 0.85],
                             binmin = 25, binmax = 60.1, binstep = 0.5)
        ax.legend()
        ax.set_xlabel("Noise [mJy beam$^{-1}$]")
        ax.set_ylabel("Count")
        ax.set_title("Continuum noise")
        pathname = os.path.join(figdir, 'dr1_cont_noise.pdf')
        plt.savefig(pathname)
        plt.close('all')

    def get_noise_table(self):
        """ 
        Get a latex formatted table of noise values
        """
        #get median and one sigma ranges for different parameters
        cont_inner = get_med_onesig(self.released['s_in_cont'])
        cont_outer = get_med_onesig(self.released['s_out_cont'])
        pol_inner = get_med_onesig(self.released['s_in_pol'])
        pol_outer = get_med_onesig(self.released['s_out_pol'])
        cube0 = get_med_onesig(self.released['rms_c0']) #in Jy
        cube1 = get_med_onesig(self.released['rms_c1'])
        cube2 = get_med_onesig(self.released['rms_c2'])

        #for polarization and line, also fidn where passed own validation
        idx_v = np.where(self.released['pass_V'] == 'True')[0]
        idx_c0 = np.where(self.released['c0_good'] == 1)[0]
        idx_c1 = np.where(self.released['c1_good'] == 1)[0]
        idx_c2 = np.where(self.released['c2_good'] == 1)[0]
        
        pol_inner_pass = get_med_onesig(self.released['s_in_pol'][idx_v])
        pol_outer_pass = get_med_onesig(self.released['s_out_pol'][idx_v])
        cube0_good = get_med_onesig(self.released['rms_c0'][idx_c0])
        cube1_good = get_med_onesig(self.released['rms_c1'][idx_c1])
        cube2_good = get_med_onesig(self.released['rms_c2'][idx_c2])
    
        #now set up the table
        #might be easiest to format this directly myself
        tablepath = os.path.join(tabledir,'dr1_noise_vals.tex')
        with open(tablepath, "w") as f:
            #table preamble
            f.write( ("\\begin{table} \n"
                      "\\centering \n"
                      "\\caption{Typical noise values} \n"
                      "\\label{tab:noise} \n"
                      "\\renewcommand{\\arraystretch}{1.2} \n"
                      "\\begin{tabular}{lll} \n"
                      "\\hline \\hline \n") )
            #table columns
            f.write( ("Data product & Released & Passed\\tablefootmark{a} \\\ \n"
                      "\\hline \n") )
            #continuum inner
            f.write( ("Continuum, inner ($\\mu$Jy bm$^{-1}$) & "
                      f"${cont_inner[0]:4.1f}^{{+{cont_inner[2]:3.1f}}}_"
                      f"{{-{cont_inner[1]:3.1f}}}$  &"
                      " -- \\\ \n ") )
            f.write( ("Continuum, outer ($\\mu$Jy bm$^{-1}$) & "
                      f"${cont_outer[0]:4.1f}^{{+{cont_outer[2]:3.1f}}}_"
                      f"{{-{cont_outer[1]:3.1f}}}$ &"
                      " -- \\\ \n ") )
            f.write( ("Stokes V, inner ($\\mu$Jy bm$^{-1}$) & "
                      f"${pol_inner[0]:4.1f}^{{+{pol_inner[2]:3.1f}}}_"
                      f"{{-{pol_inner[1]:3.1f}}}$  &"
                      f" ${pol_inner_pass[0]:4.1f}^"
                      f"{{+{pol_inner_pass[2]:3.1f}}}_"
                      f"{{-{pol_inner_pass[1]:3.1f}}}$  \\\ \n ") )
            f.write( ("Stokes V, outer ($\\mu$Jy bm$^{-1}$) & "
                      f"${pol_outer[0]:4.1f}^{{+{pol_outer[2]:3.1f}}}_"
                      f"{{-{pol_outer[1]:3.1f}}}$  &"
                      f" ${pol_outer_pass[0]:4.1f}^"
                      f"{{+{pol_outer_pass[2]:3.1f}}}_"
                      f"{{-{pol_outer_pass[1]:3.1f}}}$  \\\ \n ") )
            #put cube balues in mJY from Jy
            f.write( ("Cube0 (mJy bm$^{-1}$ ) & "
                      f"${1e3*cube0[0]:4.2f}^{{+{1e3*cube0[2]:4.2f}}}_"
                      f"{{-{1e3*cube0[1]:4.2f}}}$  &"
                      f" ${1e3*cube0_good[0]:4.2f}^"
                      f"{{+{1e3*cube0_good[2]:4.2f}}}_"
                      f"{{-{1e3*cube0_good[1]:4.2f}}}$  \\\ \n ") )
            f.write( ("Cube1 (mJy bm$^{-1}$ ) & "
                      f"${1e3*cube1[0]:4.2f}^{{+{1e3*cube1[2]:4.2f}}}_"
                      f"{{-{1e3*cube1[1]:4.2f}}}$ &"
                      f" ${1e3*cube1_good[0]:4.2f}^"
                      f"{{+{1e3*cube1_good[2]:4.2f}}}_"
                      f"{{-{1e3*cube1_good[1]:4.2f}}}$ \\\ \n ") )
            f.write( ("Cube2 (mJy bm$^{-1}$ )& "
                      f"${1e3*cube2[0]:4.2f}^{{+{1e3*cube2[2]:4.2f}}}_"
                      f"{{-{1e3*cube2[1]:4.2f}}}$ &"
                      f" ${1e3*cube2_good[0]:4.2f}^"
                      f"{{+{1e3*cube2_good[2]:4.2f}}}_"
                      f"{{-{1e3*cube2_good[1]:4.2f}}}$  \\\ \n ") )
            f.write( ("\\hline \n"
                      "\\end{tabular}\n"
                      "\\tablefoot{ \n"
                      "\\tablefootmark{a}{Good quality for cubes} } \n"
                      "\\end{table}\n" ) )
