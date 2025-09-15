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
import datetime as dt
import matplotlib.dates as mdates
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
        self.beaminfo['Field'] = np.empty(beam_table_length, dtype = "S14")
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

        # There are comments coming in with the HI table that annoy me later.
        # So get rid of them (from beam table) here:
        self.beaminfo.meta = None

    def summary(self):
        """ Print key info to screen """
        print(f'There are {len(self.obsinfo)} observations')
        print(f'That implies {40*len(self.obsinfo)} individual beams')
        print(f'The length of the continuum validation table is {len(self.continfo)}')
        print(f'The length of the polarizaiton validation table is {len(self.polinfo)}')
        print(f'The length of the cube validation table is {len(self.hiinfo)}')
        print(f'The length of the combined beam info table is {len(self.beaminfo)}')

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

        # Tom changed validation criteria slightly
        # Biggest difference is no requirement on s_in
        # Add this back here manually, so that I can match DR1
        ind_sin_fail = np.where(self.beaminfo['s_in_cont'] >= 60.)[0]
        self.beaminfo['pass'][ind_sin_fail] = 'False'
        # manually overwrite three beam ids to make things match
        ind_3 = np.where(self.beaminfo['BeamID'] == '200306072_22')[0]
        ind_2 = np.where(self.beaminfo['BeamID'] == '191004041_20')[0]
        ind_1 = np.where(self.beaminfo['BeamID'] == '190823042_11')[0]
        # 1 and 2 have sigma_in = 60 and originally passed
        self.beaminfo['pass'][ind_1] = 'True'
        self.beaminfo['pass'][ind_2] = 'True'
        # ind_3 is weird statistics, failed originally
        self.beaminfo['pass'][ind_3] = 'False'

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

    def get_cube_info(self):
        """
        Print information about released cubes
        """
        n_c2good = len(np.where(self.released['c2_good'] == 1)[0])
        n_c2okay = len(np.where(self.released['c2_ok'] == 1)[0])
        n_c1good = len(np.where(self.released['c1_good'] == 1)[0])
        n_c1okay = len(np.where(self.released['c1_ok'] == 1)[0])
        n_c0good = len(np.where(self.released['c0_good'] == 1)[0])
        n_c0okay = len(np.where(self.released['c0_ok'] == 1)[0])

        n_beams = len(self.released['c2_good'])
        n_c2bad = n_beams - n_c2good - n_c2okay
        n_c1bad = n_beams - n_c1good - n_c1okay
        n_c0bad = n_beams - n_c0good - n_c0okay

        n_good = 2*n_c2good + n_c1good + n_c0good #cube3 quality goes w/ 2
        n_okay = 2*n_c2okay + n_c1okay + n_c0okay
        n_bad = 2*n_c2bad + n_c1bad + n_c0bad
        n_bad_v2 = 4*n_beams - n_good - n_okay

        print(f"The number of good cubess is {n_good}")
        print(f"The number of okay cubes is {n_okay}")
        print(f"The number of bad cubes is {n_bad} ({n_bad_v2})")
        print(f"The total number of cubes is {4*n_beams}")

    def get_cont_csv(self):
        """
        Get csv file of DR1 released continuum beams
        """
        col_names = ['ObsID', 'Name', 'Beam', 'RA', 'Dec', 'sigma_in',
                     'sigma_out', 'bmin', 'R', 'Ex-2', 'Neg10']
        ascii.write(self.released['ObsID', 'Field', 'beam', 'RA', 'Dec', 's_in_cont', 's_out_cont',
                                  'bmin_cont', 'rat_cont', 'Ex-2_cont', 'rusc-'],
                    os.path.join(tabledir, 'dr_year1_cont.csv'),
                    format='csv',
                    overwrite=True,
                    names=col_names,
                    formats={'RA': '10.6f', 'Dec': '9.6f'}#,
                    #comment = '#'
                    )

    def get_hi_source_valid(self,
                            apertif_cat_file = os.path.join(
                                filedir, "apertif_2019_cat.v1.txt"),
                            alfalfa_cat_file = ("/Users/adams/data/alfalfa/"
                                                "a100.code12.table2.190808.csv")):
        """ 
        Validate HI source properties by undertaking a cross-matching to ALFALFA
        """
        #first specify task_ids for which comparison will be done
        #this list comes from Helga, presume it identifies observations
        #within the first source catalog that actually overlap w/ ALFALFA
        task_ids = [190914041, 190920041, 190922041, 191009039,
                    191013041, 191026001, 191027043, 191103033,
                    191115040, 191122035, 191123047, 191124035,
                    191209025, 191209026, 191219001, 191220017,
                    191222001, 191225014, 191225015, 191227013,
                    191228041]            

        apertif_cat = ascii.read(apertif_cat_file, header_start=1)

        #read in ALFALFA 100. Find this easier than constantly querying
        alfalfa100 = ascii.read(alfalfa_cat_file)

        #create empty lists to hold the information I will want
        w50_ap = []
        w50_al = []
        vhel_ap = []
        vhel_al = []
        sint_ap = []
        sint_al = []
        vali = []
        flag = []
        spat_sep = []

        #set up skycoordinates for both catalogs
        aper_coord = SkyCoord(apertif_cat['name'],unit=(u.hourangle,u.deg))
        alfalfa_coord = SkyCoord(alfalfa100['RAdeg_OC'],
                                 alfalfa100['DECdeg_OC'], unit=u.deg)
        
        #iterate through taskid list, finding sources from that taskid
        #also want to require that they are released !
        for tid in task_ids:
            ind_sources = np.where(apertif_cat['taskid'] == tid)[0]
            hi_sources = apertif_cat[ind_sources]
            #then iterate through the sources to find ALFALFA match
            for i, source in enumerate(hi_sources):
                #first check if source is released:
                beamid = f"{source['taskid']}_{source['beam']:02d}"
                if beamid in self.released['BeamID']:
                    #find the ALFALFA match
                    #take the closest spatial match as right galaxy
                    idx, sep, dist = aper_coord[
                        ind_sources][i].match_to_catalog_sky(
                            alfalfa_coord)
                    #if sep less than 1', take as match
                    if sep < 1*u.arcmin:
                        #add values to lists
                        w50_ap.append(source['w50'])
                        w50_al.append(alfalfa100['W50'][idx])
                        vhel_ap.append(source['v_sys'])
                        vhel_al.append(alfalfa100['Vhelio'][idx])
                        #need to convert flux to Jy km/s !!!
                        #channel size is 36.6 kHz
                        #this is roughly 8 km/s (but not exactly)
                        #to do this properly, use units
                        #ALFALFA is in optical, so use that
                        d_freq = 36.6 #kHz
                        d_vel = get_chan_vel(source['v_sys'], d_freq)
                        sint_ap.append(source['SJyHz'] / d_freq * d_vel)
                        print(source['v_sys'],d_vel)
                        sint_al.append(alfalfa100['HIflux'][idx])
                        spat_sep.append(sep.to(u.arcsec).value)
                        #find valid code
                        ind_release = np.where(self.released['BeamID']
                                               == beamid)[0]
                        #if source
                        flag.append(source['flag'])

        print(len(w50_ap),len(w50_al))

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
        ax.set_title("Continuum noise", size=24)
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

#helper functions
def get_chan_vel(vel, d_freq):
    rest_freq = 1420.405752*u.MHz
    hi_opt_equiv = u.doppler_optical(rest_freq)
    hi_radio_equiv = u.doppler_radio(rest_freq)
    vsys = vel*u.km/u.s
    freq_sys = vsys.to(u.Hz, equivalencies = hi_opt_equiv)
    offset_freq = freq_sys + d_freq*u.kHz
    offset_vel = offset_freq.to(u.km/u.s, equivalencies = hi_opt_equiv)
    delta_vel = offset_vel-vsys
    return np.abs(delta_vel.value)


class DR2(Beams):
    """
    Child class focused on DR2

    Attributes
    ----------

    Methods
    -------

    """

    def __init__(self,
                 obsfile=os.path.join(filedir, 'obsatdb.csv'),
                 happilifile=os.path.join(filedir, 'happili.csv'),
                 contfile=os.path.join(filedir, "cont_allbeams.csv"),
                 hifile=os.path.join(filedir, "line_allbeams.csv"),
                 polfile=os.path.join(filedir, "pol_allbeams.csv")
                 ):
        Beams.__init__(self, obsfile, happilifile, contfile, hifile, polfile)
        # All beams / observations are in DR2
        # So just need full Beams initialization
        # Obs and happili files already limit to survey observations only
        # Want to add/specify Phase1/2 here because that will generally be useful
        # Think the last phase 1 obs is 210117026 but I need to verify this!!!!
        splitid = 210117026
        self.beaminfo['Phase'] = [1 if obsid <= splitid else 2 for
                                  obsid in self.beaminfo['ObsID']]
        # add YYMM column for time-based views
        self.beaminfo['YYMM'] = [int(str(x)[:4]) for x in self.beaminfo['ObsID']]
        #self.beaminfo['YYMMDD'] = [int(str(x)[:6]) for x in self.beaminfo['ObsID']]
        self.beaminfo['date'] = [dt.date(2000+int(str(x)[:2]), int(str(x)[2:4]), int(str(x)[4:6]))
                                 for x in self.beaminfo['ObsID']]

    def get_missing_cont(self):
        """ Find and report about missing continuum validation """
        ind_no_cont = np.argwhere(np.isnan(self.beaminfo['s_in_cont']))
        # should check jointly for nans and masked values
        # but cont data has no mask
        print(f"There are {len(ind_no_cont)} beams missing continuum validation")
        print(f" (as identified by nan s_in_cont value).")
        print(f"The indices are contained in the missing_cont_beams attribute")
        self.missing_cont_beams = ind_no_cont

    def get_missing_pol(self):
        """ Find and report about missing continuum validation """
        ind_no_pol = np.argwhere(np.isnan(self.beaminfo['s_in_pol']))
        # should check jointly for nans and masked values
        # but cont data has no mask
        print(f"There are {len(ind_no_pol)} beams missing polarization validation")
        print(f" (as identified by nan s_in_pol value).")
        print(f"The indices are contained in the missing_pol_beams attribute")
        self.missing_pol_beams = ind_no_pol

    def get_missing_line(self):
        """ Find and report about missing continuum validation """
        ind_no_line = np.argwhere(np.isnan(self.beaminfo['rms_c2']))
        if len(ind_no_line) == 0:
            # check masked isnstead of na
            print('Checking for masked values')
            ind_no_line = np.where(self.beaminfo['rms_c2'].mask)[0]
        print(f"There are {len(ind_no_line)} beams missing line validation")
        print(f" (as identified by masked rms_c2 value).")
        print(f"The indices are contained in the missing_line_beams attribute")
        self.missing_line_beams = ind_no_line

    def get_cont_csv(self):
        """
        Get csv file of continuum beams
        """
        col_names = ['ObsID', 'Name', 'Beam', 'RA', 'Dec', 'Pass', 'sigma_in',
                     'sigma_out', 'bmin', 'R', 'Ex-2', 'Neg10']
        ascii.write(self.beaminfo['ObsID', 'Field', 'beam', 'RA', 'Dec', 'pass', 's_in_cont', 's_out_cont',
                                  'bmin_cont', 'rat_cont', 'Ex-2_cont', 'rusc-'],
                    os.path.join(tabledir, 'dr2_cont.csv'),
                    format='csv',
                    overwrite=True,
                    names=col_names,
                    formats={'RA': '10.6f', 'Dec': '9.6f'}#,
                    #comment = '#'
                    )

    def get_line_csv(self):
        """
        Get csv file of line valid for DR2 beams
        """
        col_names = ['ObsID', 'Name', 'Beam', 'RA', 'Dec',
                    'cube2_qual', 'cube1_qual', 'cube0_qual', 'sigma_c2',
                    'sigma_c1', 'sigma_c0',
                    'f_ex_c2', 'f_ex_c1', 'f_ex_c0',
                    'p_0.8_c2', 'p_0.8_c1', 'p_0.8_c0']

        # put sigma into mJy, rather than Jy units before formatting
        self.beaminfo['rms_c1_mJy'] = 1e3 * self.beaminfo['rms_c1']
        self.beaminfo['rms_c2_mJy'] = 1e3 * self.beaminfo['rms_c2']
        self.beaminfo['rms_c0_mJy'] = 1e3 * self.beaminfo['rms_c0']

        ascii.write(self.beaminfo['ObsID', 'Field', 'beam', 'RA', 'Dec',
                                     'c2', 'c1', 'c0', 'rms_c2_mJy',
                                     'rms_c1_mJy', 'rms_c0_mJy', 'lgfrac_c2', 'lgfrac_c1',
                                     'lgfrac_c0', 'prom_c2', 'prom_c1', 'prom_c0'],
                        os.path.join(tabledir,'dr2_line.csv'),
                        format='csv',
                        overwrite=True,
                        names=col_names,
                        formats={'sigma_c2': '4.2f', 'sigma_c1': '4.2f',
                                 'sigma_c0': '4.2f', 'f_ex_c2': '5.2f',
                                 'f_ex_c1': '5.2f', 'f_ex_c0': '5.2f',
                                 'p_0.8_c2': '4.2f', 'p_0.8_c1': '4.2f',
                                 'p_0.8_c0': '4.2f',
                                 'RA': '10.6f', 'Dec': '9.6f'}
                        )

    # make pol csv table
    def get_pol_csv(self):
        """
        csv formatted table of polarization for data release
        add ra,dec compared to paper version
        """
        col_names = ['ObsID', 'Name', 'Beam', 'RA', 'Dec', 'V valid', 'QU valid',
                     'sigma_in', 'sigma_out', 'FT_max', 'peak_inner', 'bmin',
                     'Q_beam_frac', 'U_beam_frac', 'Q_noise_frac', 'U_noise_frac']
        # col_names = ['ObsID','Name','Beam','RA','Dec','V valid','QU valid',
        #             'pol_s_in', 'pol_s_out',]
        ascii.write(self.beaminfo['ObsID', 'Field', 'beam', 'RA', 'Dec',
                                 'pass_V', 'pass_QU',
                                 's_in_pol', 's_out_pol',
                                 'ftmax', 'peak_in', 'bmin_pol', 'Q_bm_fg',
                                 'U_bm_fg', 'Q_st_fg', 'U_st_fg'],
                    os.path.join(tabledir, 'dr2_pol.csv'),
                    format='csv',
                    overwrite=True,
                    names=col_names,
                    formats={'RA': '10.6f', 'Dec': '9.6f'}
                    )

    def get_cont_valid(self):
        """
        Get the info and figures I want about cont valid for DR2
        """
        # First, want to get median noise values for both phases, both all and valid only
        ind1 = np.where(self.beaminfo['Phase'] == 1)[0]
        ind2 = np.where(self.beaminfo['Phase'] == 2)[0]
        beam1 = self.beaminfo[ind1]
        beam2 = self.beaminfo[ind2]
        valid1 = np.where(beam1['pass'] == 'True')[0]
        valid2 = np.where(beam2['pass'] == 'True')[0]
        beam1_valid = beam1[valid1]
        beam2_valid = beam2[valid2]
        inner1 = np.nanmedian(beam1['s_in_cont'])
        inner1_valid = np.nanmedian(beam1_valid['s_in_cont'])
        outer1 = np.nanmedian(beam1['s_out_cont'])
        outer1_valid = np.nanmedian(beam1_valid['s_out_cont'])
        inner2 = np.nanmedian(beam2['s_in_cont'])
        inner2_valid = np.nanmedian(beam2_valid['s_in_cont'])
        outer2 = np.nanmedian(beam2['s_out_cont'])
        outer2_valid = np.nanmedian(beam2_valid['s_out_cont'])
        print(f"Median inner noise for Phase 1 valid"
              f" (alL) cont is {inner1_valid:4.1f} ({inner1:4.1f})")
        print(f"Median outer noise for Phase 1 valid"
              f" (alL) cont is {outer1_valid:4.1f} ({outer1:4.1f})")
        print(f"Median inner noise for Phase 2 valid"
              f" (alL) cont is {inner2_valid:4.1f} ({inner2:4.1f})")
        print(f"Median outer noise for Phase 2 valid"
              f" (alL) cont is {outer2_valid:4.1f} ({outer2:4.1f})")
        # Also print an overview of the number of beams
        print(f"There are {len(beam1)} possible beams in Phase 1, {len(beam1_valid)} of which pass validation")
        print(f"There are {len(beam2)} possible beams in Phase 2, {len(beam2_valid)} of which pass validation")
        # But here I need to account for the fact that not all beams exist (i.e., processing failed)
        nodata1 = np.argwhere(np.isnan(beam1['s_in_cont']))
        nodata2 = np.argwhere(np.isnan(beam2['s_in_cont']))
        print(f"There are {len(beam1)-len(nodata1)} released beams in Phase 1, {len(beam1_valid)} of which pass validation")
        print(f"There are {len(beam2) - len(nodata2)} released beams in Phase 2, {len(beam2_valid)} of which pass validation")
        # Plot histograms of the inner noise
        fig, ax = plot_hist(beam1['s_in_cont'], beam2['s_in_cont'],
                            beam1_valid['s_in_cont'], beam2_valid['s_in_cont'],
                            colors=['#AA3377', '#4477AA','#AA3377', '#4477AA'],
                            alpha = [0.5, 0.5, 1, 1],
                            labels=['Phase 1 all', 'Phase 2 all','Phase 1 valid','Phase 2 valid'],
                            #alpha=[1.0, 0.85],
                            binmin=20, binmax=100.1, binstep=0.5)
        ax.legend()
        ax.set_xlabel("Noise [mJy beam$^{-1}$]")
        ax.set_ylabel("Count")
        ax.set_title("Inner continuum noise", size=24)
        pathname = os.path.join(figdir, 'dr2_cont_inner_noise.pdf')
        plt.savefig(pathname)
        plt.close('all')
        # Plot histograms of the outer noise
        fig, ax = plot_hist(beam1['s_out_cont'], beam2['s_out_cont'],
                            beam1_valid['s_out_cont'], beam2_valid['s_out_cont'],
                            colors=['#AA3377', '#4477AA','#AA3377', '#4477AA'],
                            alpha = [0.5, 0.5, 1, 1],
                            labels=['Phase 1 all', 'Phase 2 all','Phase 1 valid','Phase 2 valid'],
                            #alpha=[1.0, 0.85],
                            binmin=20, binmax=100.1, binstep=0.5)
        ax.legend()
        ax.set_xlabel("Noise [mJy beam$^{-1}$]")
        ax.set_ylabel("Count")
        ax.set_title("Outer continuum noise", size=24)
        pathname = os.path.join(figdir, 'dr2_cont_outer_noise.pdf')
        plt.savefig(pathname)
        plt.close('all')
        # get time-based views of what is happening
        # want to sort things based on YYMM column
        yymm_array = np.unique(self.beaminfo['YYMM'])
        s_in_cont_yymm_median = [np.nanmedian(self.beaminfo['s_in_cont'][self.beaminfo['YYMM']==t]) for t in yymm_array]
        yymm_date_array = [dt.date(2000+int(str(t)[0:2]), int(str(t)[2:4]), 15) for t in yymm_array]

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes((0.11, 0.1, 0.85, 0.85))
        ax.scatter(self.beaminfo['date'],self.beaminfo['s_in_cont'], marker='.', color='gray')
        ax.scatter(yymm_date_array, s_in_cont_yymm_median, s=50, marker='s', color='red')
        # Major ticks every half year, minor ticks every month,
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
        ax.xaxis.set_minor_locator(mdates.MonthLocator())
        ax.set_xlabel('Date')
        ax.set_ylabel('Inner noise [mJy beam$^{-1}$]')
        ax.set_ylim(15,250)
        pathname=os.path.join(figdir, 'dr2_cont_inner_noise_time.pdf')
        plt.savefig(pathname)
        plt.close('all')

