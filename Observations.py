#Observation object

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Object to hold Apertif observational info
Children objects for more specific purposes:
DR1 : First data release paper
Census: Census (early Nov 2021) for finalizing to end of ops
"""

import os
import matplotlib.pyplot as plt
from astropy.io import ascii
import numpy as np
from astropy.table import Table, join
from functions.plots import plot_sky_view


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


class Observations(object):
    """
    A class holding Apertif observational information

    ...

    Attributes
    ----------
    obsinfo : Astropy Table
        table of observational information

    Methods
    -------
    summary

    """
    def __init__(self,
                 obsfile = os.path.join(filedir,'obsatdb.csv')):
        
        """
        Construct initial attributes for Observations object

        Parameters
        ----------
        obsfile : str
            Full path to atdb observational file
        """
        self.obsinfo = ascii.read(obsfile)

    def summary(self):
        """
        Print a summary of observations
        """
        #get number of medium-deep fields
        ind_ames = [i for i, s in enumerate(self.obsinfo['name']) if 'M' in s]
        ames = self.obsinfo[ind_ames]
        #get unique & count
        unique_ames, count_ames = np.unique(ames['name'], return_counts = True)
        #limit to those with more than two fields
        n_rep_ames = 0
        n_mult_obs_ames = 0
        for field, count in zip(unique_ames, count_ames):
            if count > 2:
                print( ("There are {0} observations of "
                        "field {1}").format(count,field) )
                n_rep_ames = n_rep_ames + 1
                n_mult_obs_ames = n_mult_obs_ames + count

        #get number of shallow fields
        ind_awes = [i for i, s in enumerate(self.obsinfo['name']) if 'S' in s]
        awes = self.obsinfo[ind_awes]
        #get number unique
        unique_awes, count_awes = np.unique(awes['name'], return_counts = True)
        #get number w/ repeats
        n_mult_visit = len(np.where(count_awes > 1)[0])

        print(("There are {0} observations of "
               "{1} independent medium-deep fields").format(len(ames),
                                                            len(unique_ames)))

        print(("There are {0} medium-deep fields "
               "with repeated coverage, over a total of "
               "{1}  observations").format(n_rep_ames,
                                           n_mult_obs_ames))

        print(("There are {0} observations of "
               "{1} independent wide/shallow fields").format(len(awes),
                                                             len(unique_awes)))
        print(("There are {0} wide fields "
               "with repeat observations").format(n_mult_visit))

        print(("There are {0} observations in total").format(len(self.obsinfo)))

        print( ("There are {0} unique fields in total, "
                "for a sky coverage of {1:4.0f} square degrees").format(
                    len(np.unique(self.obsinfo['name'])),
                    6.44*len(np.unique(self.obsinfo['name']))) )

 

class Census(Observations):
    """
    Child class focused on census of good/bad observations

    ...

    Attributes
    ----------
    obsinfo : Astropy Table
        table of observational information
    censusinfo : Astropy Table
        table of census quality

    Methods
    -------
    plot_obs_census
    get_census_details
    print_census_summary
    """
    def __init__(self,
                 obsfile = os.path.join(filedir,'obsatdb.csv'),
                 censusfile = os.path.join(filedir,'obscensus.csv')):
        """
        Construct initial attributes for Census object

        Parameters
        ----------
        obsfile : str
            Full path to atdb observational file
        censusfile : str
            Full path to census file
        """
        Observations.__init__(self,obsfile)
        self.censusinfo = ascii.read(censusfile, format='csv',
                                     comment="\s*#",)

        #add frequency information to census
        #do this based on data
        self.censusinfo['Freq_cen'] = np.where(
            ( self.censusinfo['taskID'] > 210118000 ), 1370.0, 1280.0 )

        #now join the observation to census info
        #everything that is specific to Census class
        #will be done w/ census info
        self.censusinfo = join(self.obsinfo,self.censusinfo,
                               keys = 'taskID', join_type='outer')

        #get field-based census
        #get unique fields
        fields, field_inds = np.unique( self.censusinfo['name_1'],
                                        return_index = True)

        #set up arrays to hold other info
        telescopes = np.empty(len(fields),dtype=object) #will hold telescope arrays
        ndishes = np.empty(len(fields)) # will be number of dishes
        nobs = np.empty(len(fields)) #number of obs
        avg_dishes = np.empty(len(fields)) #number of dishes / obs
        mid_freq = np.empty(len(fields)) #center frequency of all obs (mean)
        frac_before_1oct19 = np.empty(len(fields)) #fraction of obs before 1oct2019
        frac_before_11dec19 = np.empty(len(fields)) #fraction before 11dec19
        frac_cd = np.empty(len(fields)) #fraction of dishes that are CD per obs
        N_CD = np.empty(len(fields))

        #iterate through each field
        for i,f in enumerate(fields):
            f_inds = np.where( self.censusinfo['name_1'] == f )[0]
            telescope_list = [list(x) for x in
                              self.censusinfo['goodtelescopes'][f_inds] ]
            #add to telescopes array as a list
            telescopes[i] = [num for elem in telescope_list for num in elem]
            ndishes[i] = len(telescopes[i])
            nobs[i] = len(f_inds)
            avg_dishes[i] = ndishes[i] / nobs[i]
            mid_freq[i] = np.mean(self.censusinfo['Freq_cen'][f_inds])
            N_C = telescopes[i].count('C')
            N_D = telescopes[i].count('D')
            N_CD[i] = N_C + N_D
            frac_cd[i] = N_CD[i] / nobs[i]

        #make a table and add it as an attribute
        self.field_census = Table()
        self.field_census['Field'] = fields
        self.field_census['RA'] = self.censusinfo['field_ra'][field_inds]
        self.field_census['Dec'] = self.censusinfo['field_dec'][field_inds]
        self.field_census['telescopes'] = telescopes
        self.field_census['Ndishes'] = ndishes
        self.field_census['Nobs'] = nobs
        self.field_census['AvgDishes'] = avg_dishes
        self.field_census['Mid_freq'] = mid_freq
        self.field_census['Avg_CD_obs'] = frac_cd
        self.field_census['N_CD'] = N_CD

        

        
    def plot_obs_census(self, view,
                        surveypointings = os.path.join(
                            filedir,
                            'all_pointings.v7.18jun20.txt') ):
        """
        Plot views of observational census

        Want to plot with a colorbar scale
        And provide the view as a keyword argument

        Parameters
        ----------
        view : str
           The view to be plotted. Options are:
           ["Nobs", "avg_cd_obs", "avg_dishes"]
        surveypointings : str
           Full path to file of survey pointings to be plotted for comparison
        """
        if view not in ["Nobs", "avg_cd_obs", "avg_dishes"]:
            print("View name not found.")
            print("Defaulting to number of observations")
            view = "Nobs"
        
        if view == "Nobs":
            ind_1 = np.where(self.field_census['Nobs'] == 1)[0]
            ind_2 = np.where(self.field_census['Nobs'] == 2)[0]
            ind_3 = np.where(self.field_census['Nobs'] == 3)[0]
            ind_10 = np.where(self.field_census['Nobs'] >= 10)[0]
            ind_4_9 = np.where( np.logical_and(
                self.field_census['Nobs'] >= 4,
                self.field_census['Nobs'] <= 9) )[0]
            ralist = [self.field_census['RA'][ind_1],
                      self.field_census['RA'][ind_2],
                      self.field_census['RA'][ind_3],
                      self.field_census['RA'][ind_4_9],
                      self.field_census['RA'][ind_10]]
            declist = [self.field_census['Dec'][ind_1],
                       self.field_census['Dec'][ind_2],
                       self.field_census['Dec'][ind_3],
                       self.field_census['Dec'][ind_4_9],
                       self.field_census['Dec'][ind_10] ]
            labellist = ["1 visit", "2 visits", "3 visits",
                         "4-9 visits", "10 or more visits"]
        if view == "avg_cd_obs":
            ind_0 = np.where(self.field_census['Avg_CD_obs'] == 0)[0]
            ind_half =  np.where( np.logical_and(
                self.field_census['Avg_CD_obs'] > 0,
                self.field_census['Avg_CD_obs'] <= 0.5 ) )[0]
            ind_1 = np.where( np.logical_and(
                self.field_census['Avg_CD_obs'] > 0.5,
                self.field_census['Avg_CD_obs'] <1 ) )[0]
            ind_2 = np.where( self.field_census['Avg_CD_obs'] >= 1)[0]
            ralist = [self.field_census['RA'][ind_0],
                      self.field_census['RA'][ind_half],
                      self.field_census['RA'][ind_1],
                      self.field_census['RA'][ind_2]]
            declist = [self.field_census['Dec'][ind_0],
                       self.field_census['Dec'][ind_half],
                       self.field_census['Dec'][ind_1],
                       self.field_census['Dec'][ind_2]]
            labellist = ['No C/D', "Avg lte 0.5 per obs",
                         "Avg lt 1 per obs", "Avg gte 1 per obs" ]

        if view == "avg_dishes":
            ind_12 = np.where(self.field_census['AvgDishes'] == 12)[0]
            ind_10_12 = np.where(np.logical_and(
                self.field_census['AvgDishes'] < 12,
                self.field_census['AvgDishes'] >=10) )[0]
            ind_8_10 = np.where(np.logical_and(
                self.field_census['AvgDishes'] < 10,
                self.field_census['AvgDishes'] > 8))[0]
            ind_8 = np.where(self.field_census['AvgDishes'] <= 8)[0]
            ralist = [self.field_census['RA'][ind_12],
                      self.field_census['RA'][ind_10_12],
                      self.field_census['RA'][ind_8_10],
                      self.field_census['RA'][ind_8]
            ]
            declist = [self.field_census['Dec'][ind_12],
                       self.field_census['Dec'][ind_10_12],
                       self.field_census['Dec'][ind_8_10],
                       self.field_census['Dec'][ind_8]
            ]
            labellist = ["Avg 12 dishes", "Avg between 10 and 12",
                         "Avg between 8-10", "Avg le 8"]


        plot_sky_view(ralist, declist, labellist, view,
                      surveypointings = surveypointings)

        

    def print_census_summary(self):
        """
        Print relevant info from census
        """
        #identify fields missing C + D
        ind_0 = np.where(self.field_census['Avg_CD_obs'] == 0)[0]
        ind_half =  np.where( np.logical_and(
            self.field_census['Avg_CD_obs'] > 0,
            self.field_census['Avg_CD_obs'] <= 0.5 ) )[0]
        ind_one = np.where( np.logical_and(
            self.field_census['Avg_CD_obs'] > 0.5,
            self.field_census['Avg_CD_obs'] < 1
            ))[0]
        print( ("There are {} observations with no C or D").
               format(len(ind_0)))
        for i in ind_0:
            print("Field {} has no C or D coverage".format(
                self.field_census['Field'][i]))
        for i in ind_half:
            print("Field {} has <= 1/2 C/D on average".format(
                self.field_census['Field'][i]))
        for i in ind_one:
            print(("Field {0} has {1} C/D on average").format(
                self.field_census['Field'][i],
                self.field_census['Avg_CD_obs'][i]))

            
    def get_fields_reobserve(self,outfile):
        """
        Get fields to reobserve

        Parameters
        ----------
        outfile : str
            Full path of file to write results to
        """
        #first find fields missing CD
        ind_nocd = np.where(self.field_census['Avg_CD_obs'] == 0)[0]
        fields_nocd = self.field_census[ind_nocd]
        fields_nocd['NoCD'] = np.full(len(fields_nocd),'True')

        #then find fields with one obs and <= 9 dishes
        ind_lackdishes = np.where( np.logical_and(
            self.field_census['Nobs'] ==1,
            self.field_census['AvgDishes'] <=9
            ))[0]
        fields_lackdishes = self.field_census[ind_lackdishes]
        fields_lackdishes['LackDishes'] = np.full(len(fields_lackdishes),
                                                  'True')

        #do an outer join of these tables
        fields_reobserve =  join(fields_nocd, fields_lackdishes,
                                 keys = 'Field', join_type='outer')

        #and save to a file
        fields_reobserve['Field','NoCD','LackDishes'].write(
            outfile, format='ascii', overwrite=True)
