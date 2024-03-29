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
import casacore.tables as pt
import glob

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
aperinfodir = this_dir[:-4]
filedir = os.path.join(aperinfodir,"files")
#print(filedir)


def check_cd(tid,Nbeams):
    """
    Take a taskid and number of beams to fail
    Check how many beams are missing C&D
    Return False if too many are missing
    Inputs:
    - tid: Taskid to check
    - Nbeams: Number of beams that equals a failure
    Outputs:
    - True/False: True, observation is good.
                  False, observation fails
    """
    #first set an array to hold output
    #for each beam
    #Set to be True, update to False if beam fails
    beam_array = np.full(40,True)
    #iterate through each beam
    for bm in range(40):
        #check that bandpass table exists
        #first get right happili data path
        datapath = get_data_path(bm)
        #then the raw directory path
        rawbeamdir = os.path.join(datapath,tid,
                                  '{0:02d}/raw'.format(bm))
        #and find bandpass, if it exists
        bptest = glob.glob(os.path.join(rawbeamdir,"*Bscan"))
        if len(bptest) == 1:
            bppath = bptest[0]
        if len(bpest) == 0:
            print('No bandpass file; marking as failed')
            beam_array[bm] = False
        else:
            print('Help! Found {} bandpass files'.format(len(bptest)))
        
        

def get_data_path(bm):
    """
    Given a beam number, return correct
    /data?/apertif, presuming on happili-01
    """
    #just in case passed as a string
    bm = int(bm)
    if bm < 10:
        datapath = '/data/apertif/'
    elif bm < 20:
        datapath = '/data2/apertif/'
    elif bm < 30:
        datapath = '/data3/apertif/'
    else:
        datapath = '/data4/apertif/'
        
    #return the path
    return datapath
    

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


def get_bad_pol():
    """
    Have to find observation/beam combinations where
    there is no pol calibrator
    Pipeline errors mean polarization data products are 
    produced,  but they're not calibrated.
    Need to find and remove from ALTA
    Also happili ?
    Two cal tables: .Kcross and .Xf indicate pol cal is present
    """
    #first, get directory list
    taskdirlist = get_dir_list()
    #then go through each, checking each beam
    #create empty lists to hold things
    beamlist = []
    tasklist = []
    for taskdir in taskdirlist:
        taskid = taskdir[-9:]
        day = np.float(taskid[0:6])
        #ignore July and start August, to make my life easier
        if day > 190806:
            for i in range(40):
                #!!!Have to updatee for right happili node
                if i < 10:
                    rawbeamdir = os.path.join(taskdir,
                                              '{0:02d}/raw'.format(i))
                elif i < 20:
                    happilidir = os.path.join('/data2/apertif/',taskid)
                    rawbeamdir = os.path.join(happilidir,
                                              '{0:02d}/raw'.format(i))
                elif i < 30:
                    happilidir = os.path.join('/data3/apertif/',taskid)
                    rawbeamdir = os.path.join(happilidir,
                                              '{0:02d}/raw'.format(i))
                else:
                    happilidir = os.path.join('/data4/apertif/',taskid)
                    rawbeamdir = os.path.join(happilidir,
                                              '{0:02d}/raw'.format(i))
                #first check dir exists, otherwise can skip completely
                if os.path.exists(rawbeamdir):
                    #need to check fluxcal also
                    #think no fluxcal should do things properly
                    #and also no pol solutions?
                    #but first verify that
                    polXf = glob.glob(os.path.join(rawbeamdir,'3C???.Xf'))
                    if len(polXf) == 0:
                        print(('No Xf pol solution for'
                               ' {0}, beam {1}').format(taskid,i))
                        beamlist.append(i)
                        tasklist.append(taskid)
    #write everything out to a file
    #make a table
    badpol = Table([tasklist,beamlist],names=['taskid','beam'])
    ascii.write(badpol,os.path.join(filedir,'badpol.csv'),format='csv',
                overwrite=True)

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
    t['taskid'] = np.full(len(tasklist),'YYMMDDTTT')
    t['fluxcal'] = np.full(len(tasklist),'3C???')
    t['fluxcal_firsttaskid'] = np.full(len(tasklist),'YYMMDDTTT')
    t['fluxcal_lasttaskid'] = np.full(len(tasklist),'YYMMDDTTT')
    t['polcal'] = np.full(len(tasklist),'3C???')
    t['polcal_firsttaskid'] = np.full(len(tasklist),'YYMMDDTTT')
    t['polcal_lasttaskid'] = np.full(len(tasklist),'YYMMDDTTT')
    t['apercal_version'] = np.full(len(tasklist),'v?.?-NNN-12345678')
    """
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
    """    

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
    ascii.write(t,os.path.join(filedir,'happili.csv'),format='csv',
                overwrite=True)



def get_cal_info(taskdir):
    """
    Helper script to retrieve calibrator information
    Inputs:
        taskdir (str): Full path for taskid on happili-01
    """
    #Could use OSA reports, but don't always exist
    #better to get directly
    #first cal is on happili-01, last on happili-04
    taskid = taskdir[-9:]
    path1 = os.path.join(taskdir,"qa/{0}_obs.ecsv".format(taskid))
    path4 = path1[0:5] + "4" + path1[5:]
    #add a try except in file doesn't exist
    #this would indicate a processing / autocal issue
    try:
        info1 = ascii.read(path1)
        info4 = ascii.read(path4)
        fluxcal = info1['Flux_Calibrator'][0]
        polcal = info1['Pol_Calibrator'][0]
    except:
        fluxcal = None
        polcal = None
    #not all files have id info, so try except
    #define variables first, w/ placeholder
    fluxid_start = None
    fluxid_end = None
    polid_start = None
    polid_end = None
    try:
        fluxstring1 = info1['Flux_Calibrator_Obs_IDs'][0]
        #Note that lists of cals may not be sorted, so do this manually
        #use helper function since I have to do four times
        #twice for first, twice for last
        fluxid_start = get_id(fluxstring1,place='first')
        fluxstring4 = info4['Flux_Calibrator_Obs_IDs'][0]
        fluxid_end = get_id(fluxstring4,place='last')
        polstring1 = info1['Pol_Calibrator_Obs_IDs'][0]
        polid_start = get_id(polstring1,place='first')
        polstring4 = info4['Pol_Calibrator_Obs_IDs'][0]
        polid_end = get_id(polstring4,place='last')
    except:
        fluxid_start = None
        fluxid_end = None
        polid_start = None
        polid_end = None

    return fluxcal, fluxid_start, fluxid_end, polcal, polid_start, polid_end

def get_id(idstring,place='first'):
    """
    Get id from the string listing all ten ids
    place = 'first' or 'last'
    """
    #get rid of spaces and split on commas   
    ids = idstring.replace(" ","").split(',')
    #sort
    ids.sort()
    #get right id
    if place == 'first':
        finalid = ids[0]
    elif place == 'last':
        finalid = ids[-1]
    else:
        print("'place' not recognized. Must be 'first' or 'last'")

    return finalid
                                            

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
    

    
