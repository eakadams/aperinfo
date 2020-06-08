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

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
aperinfodir = this_dir[:-4]
filedir = os.path.join(aperinfodir,"files")
#print(filedir)


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
    t['fluxcal'] = np.empty(len(tasklist),dtype=str)
    t['fluxcal_firsttaskid'] = np.empty(len(tasklist),dtype=str)
    t['fluxcal_lasttaskid'] = np.empty(len(tasklist),dtype=str)
    t['polcal'] = np.empty(len(tasklist),dtype=str)
    t['polcal_firsttaskid'] = np.empty(len(tasklist),dtype=str)
    t['polcal_lasttaskid'] = np.empty(len(tasklist),dtype=str)
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
        #run_times_01 = get_run_times(task)
        #skip run times for now - not critical
        #informational and want to do
        #but calibrator / apercal information is what is needed now
        (fluxcal, fluxid_start, fluxid_end,
         polcal, polid_start, polid_end) = get_cal_info(task)
        t['fluxcal'][i] = fluxcal
        t['fluxcal_firsttaskid'][i] = fluxid_start
        t['fluxcal_lasttaskid'][i] = fluxid_end
        t['polcal'][i] = polcal
        t['polcal_firsttaskid'][i] = polid_start
        t['polcal_lasttaskid'][i] = polid_end

    #write table out
    ascii.write(t,os.path.join(filedir,'happili.csv',format='csv'))



def get_cal_info(taskdir):
    """
    Helper script to retrieve calibrator information
    Inputs:
        taskdir (str): Full path for taskid on happili-01
    """
    #Could use OSA reports, but don't always exist
    #better to get directly
    #first cal is on happili-01, last on happili-04
    print(taskdir)
    taskid = taskdir[-9:]
    path1 = os.path.join(taskdir,"qa/{0}_obs.ecsv".format(taskid))
    path4 = path1[0:5] + "4" + path1[5:]
    info1 = ascii.read(path1)
    info4 = ascii.read(path4)
    print(info1.colnames)
    fluxcal = info1['Flux_Calibrator'][0]
    fluxstring1 = info1['Flux_Calibrator_Obs_IDs'][0]
    fluxid_start = fluxstring1[0:9]
    fluxstring4 = info4['Flux_Calibrator_Obs_IDs'][0]
    fluxid_end = fluxstring[-9:]
    polcal = info1['Pol_Calibrator'][0]
    polstring1 = info1['Pol_Calibrator_Obs_IDs'][0]
    polid_start = polstring1[0:9]
    polstring4 = info4['Pol_Calibrator_Obs_IDs'][0]
    polid_end = polstring[-9:]

    return fluxcal, fluxid_start, fluxid_end, polcal, polid_start, polid_end

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
    path1 = os.path.join(taskdir,'qa/apercal_performance/apercal_log_timeinfo_happili-01.csv')
    path2 = path1[0:5] + "2" + path1[5:-5] + "2" + path1[-4:]
    path3 = path1[0:5] + "2" + path1[5:-5] + "2" + path1[-4:]
    path4 = path1[0:5] + "2" + path1[5:-5] + "2" + path1[-4:]
    runtimes01 = ascii.read(path1)
    
    

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
    taskid = taskdir[-9:]
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
    

    
