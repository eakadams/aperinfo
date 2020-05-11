#Functions for getting info from happilis

from __future__ import print_function

__author__="E.A.K. Adams"

"""
Functions related to pulling information
from happili nodes.
Assumes you are running from happili-01;
the only way to see all happili nodes
"""

import glob
import sys
from astropy.table import Table
import numpy as np
from astropy.io import ascii
from datetime import datetime
import os

#should be updated for local location
#need to figure out to make this global
#and a level above
#but I'm lazy
filedir = "/home/adams/aperinfo/files/"


def get_dir_list():
    """
    Get list of taskid directories on happili-01
    This shoudl be all survey observations 
    processed to date.
    And only survey observations
    But if others sneak in, that should be okay 
    """
    taskdirlist = glob.glob('/data/apertif/?????????')
    #update, better option, from Robert
    taskdirlist = glob.glob(
    "/data/apertif/[1-2][0-9][0-1][0-9][0-3][0-9][0-9][0-9][0-9]")
    #and return it in sorted order
    taskdirlist.sort()
    return taskdirlist

def make_happili_obs_table():
    """
    Make an observation / field - based table
    of information that is on happili.
    That is, this is a table with a row per each taskid
    """
    tasklist = get_dir_list()
    #make a table that is the same length
    #create columns that want in the end
    t=Table()
    t['taskid'] = np.empty(len(tasklist),dtype=str)
    t['apercal_version'] = np.empty(len(tasklist),dtype=str)
    t['apercal_01_runtime']=np.empty(len(tasklist),dtype=str)
    t['apercal_02_runtime']=np.empty(len(tasklist),dtype=str)
    t['apercal_03_runtime']=np.empty(len(tasklist),dtype=str)
    t['apercal_04_runtime']=np.empty(len(tasklist),dtype=str)
    #could add information for each step
    #this might be useful for statistics and identifying
    #possible problems
    #so I'll go ahead and do it
    t['apercal_01_runtime_prepare'] = np.empty(len(tasklist),dtype=str)
    t['apercal_01_runtime_split'] = np.empty(len(tasklist),dtype=str)
    t['apercal_01_runtime_preflag'] = np.empty(len(tasklist),dtype=str)
    t['apercal_01_runtime_crosscal'] = np.empty(len(tasklist),dtype=str)
    t['apercal_01_runtime_convert'] = np.empty(len(tasklist),dtype=str)
    t['apercal_01_runtime_selfcal+continuum+polarisation'] = np.empty(len(tasklist),dtype=str)
    t['apercal_01_runtime_line'] = np.empty(len(tasklist),dtype=str)
    t['apercal_01_runtime_tranfer'] = np.empty(len(tasklist),dtype=str)
    t['apercal_02_runtime_prepare'] = np.empty(len(tasklist),dtype=str)
    t['apercal_02_runtime_split'] = np.empty(len(tasklist),dtype=str)
    t['apercal_02_runtime_preflag'] = np.empty(len(tasklist),dtype=str)
    t['apercal_02_runtime_crosscal'] = np.empty(len(tasklist),dtype=str)
    t['apercal_02_runtime_convert'] = np.empty(len(tasklist),dtype=str)
    t['apercal_02_runtime_selfcal+continuum+polarisation'] = np.empty(len(tasklist),dtype=str)
    t['apercal_02_runtime_line'] = np.empty(len(tasklist),dtype=str)
    t['apercal_02_runtime_tranfer'] = np.empty(len(tasklist),dtype=str)
    t['apercal_03_runtime_prepare'] = np.empty(len(tasklist),dtype=str)
    t['apercal_03_runtime_split'] = np.empty(len(tasklist),dtype=str)
    t['apercal_03_runtime_preflag'] = np.empty(len(tasklist),dtype=str)
    t['apercal_03_runtime_crosscal'] = np.empty(len(tasklist),dtype=str)
    t['apercal_03_runtime_convert'] = np.empty(len(tasklist),dtype=str)
    t['apercal_03_runtime_selfcal+continuum+polarisation'] = np.empty(len(tasklist),dtype=str)
    t['apercal_03_runtime_line'] = np.empty(len(tasklist),dtype=str)
    t['apercal_03_runtime_tranfer'] = np.empty(len(tasklist),dtype=str)
    t['apercal_04_runtime_prepare'] = np.empty(len(tasklist),dtype=str)
    t['apercal_04_runtime_split'] = np.empty(len(tasklist),dtype=str)
    t['apercal_04_runtime_preflag'] = np.empty(len(tasklist),dtype=str)
    t['apercal_04_runtime_crosscal'] = np.empty(len(tasklist),dtype=str)
    t['apercal_04_runtime_convert'] = np.empty(len(tasklist),dtype=str)
    t['apercal_04_runtime_selfcal+continuum+polarisation'] = np.empty(len(tasklist),dtype=str)
    t['apercal_04_runtime_line'] = np.empty(len(tasklist),dtype=str)
    t['apercal_04_runtime_tranfer'] = np.empty(len(tasklist),dtype=str)
    
    #once all columns are defined
    #iterate through and fill the table
    for i,task in enumerate(tasklist):
        #get task id
        tid = task[-9:]
        t['taskid'][i] = tid
        #get apercal version
        #only check happili-01 here
        apercal_vers = get_apercal_version(task)
        t['apercal_version'][i]=apercal_vers
        #get apercal run times
        #do to read csv on each happili, then have all info
        run_times_01 = get_run_times(task)

def get_run_times(taskdir):
    """
    Helper script to retrieve run times
    Actually - let's be smart and have this script
    retrieve from all nodes
    and return a 2-d array (or table?) that has
    information for each run time type, plus happili node
    Inputs:
        taskdir (str): Full path for taskid (on happili-01)
    Outputs:
        runtimes (nparrary): 4x9 array with run time information
    """
    
    

def get_apercal_version(taskdir):
    """
    Based on code by R. Schulz. 
    Takes a taskid full directory path
    And return apercal version string
    Inputs:
         taskdir (str): Full directory path to a taskid
    Outputs:
         apercal_vers (str): Apercal version. None if no version
    """
    logfile_name = os.path.join(taskdir, "apercal.log")
    search_phrase = "Apercal version:"
    # make sure the logfile exists
    if not os.path.exists(logfile_name):
        print("WARNING: Did not find logfile for {}".format(
            os.path.basename(taskdir)))
        apercal_vers = None
    else:
        # read logfile
        with open(logfile_name, "r") as logfile:
            # go through the lines
            for logline in logfile:

                # check for interesting line
                if search_phrase in logline:
                    apercal_vers = logline.split(": ")[-1].split("\n")[0]
                    print("{0} {1}".format(os.path.basename(taskid),
                                           apercal_vers))
                                           #logline.split(": ")[-1].split("\n")[0]))

    return apercal_vers
    
def make_happili_beam_table():
    """
    Make a beam-based table based on happili information
    This is a longer table that has an entry per taskid/beam pair
    """
    
def update_happili_obstable(oldobsfile):
    """
    Take an older observation file
    Only update for the new taskids
    Need to test/confirm this is faster 
    than regenerating everything
    """
    
def update_happili_beamtable(oldbeamfile):
    """
    Take an older beam-based file
    Only update for the new taskids
    Need to test/confirm this is faster 
    than regenerating everything
    """
    
