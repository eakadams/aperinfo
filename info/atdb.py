#Functions for getting info from ATDB

from __future__ import print_function

__author__="E.A.K. Adams"

"""
Functions related to querying ATDB for
retrieving observational information.
Relies on atdbquery package from
V.A. Moss (cosmicpudding). Should be 
in python path or added explicitly.
Currently using a modified version
with:
valid_statuses = "removed"
This makes the query significantly faster
"""

import sys
sys.path.append('/home/adams/atdbquery')
import atdbquery
from astropy.table import Table
import numpy as np
from astropy.io import ascii
from datetime import datetime
import os

#should be updated for local location
filedir = "/home/adams/aperinfo/files/"

def get_obslist():
    """
    Query ATDB to get observation list
    Simple but slow
    If I use this code more than once in module,
    I'm doing something wrong
    """
    obslist = atdbquery.atdbquery('imaging',False,False)
    return obslist

def get_obstable(write=True):
    """
    Get a table object that has basic observational information I want
    """
    #get info from ATDB
    #don't really need separate code from this
    #should be only time I call it
    obslist = get_obslist()
    #turn into a table; this is fast
    fulltable = Table(obslist)
    #then select based on duration for field scans only
    #take longer than 10 hours
    ind_long = np.where(fulltable['duration'] >= 36000.)[0]
    targettable = fulltable[ind_long]
    #that is table of all long (presumably target) observations
    #but could include some test observations or early science
    #want to sort on name to find those that are survey fields
    #know there is at least one exception
    #M0155+3622 191030203 was a test with calculated beamweights
    #iterate through because I'm lazy
    #i also want to first sort by taskid because I want things
    #to come in the right order
    targettable.sort('taskID')
    surveyinds = []
    otherinds = []
    testinds = []
    earlyscienceinds = []
    argoinds = []
    failedinds = []
    for i,(name,taskid,quality) in enumerate(targettable['name','taskID','quality']):
            #first grab specific fields that I know are tests
            if ( (taskid == '191030203') or ('test' in name) or
                 (name == 'L1035+5720') or (taskid == '191216159') ):
                #known test observations
                testinds.append(i)
            elif name[0:4] == 'ARGO':
                #find ARGO fields
                argoinds.append(i)
            elif (int(taskid) < 190702000) and (int(taskid) > 190409000):
                #find early science fields based on date
                earlyscienceinds.append(i)
            elif ( taskid == '210718041' ):
                #observation that accidentally deleted
                failedinds.append(i)
            elif ((len(name) == 10) and (int(taskid) > 190702000)
                  and ((quality == 'good') or (quality=='unknown')) ):
                #find survey fields based on length and date
                surveyinds.append(i)
            elif ((len(name) == 10) and (int(taskid) > 190702000)
                  and ((quality == 'bad')) ):
                #those that are marked bad
                failedinds.append(i)
            else:
                #dump everything else somewhere
                otherinds.append(i)

    #get various table subsets
    surveyfields = targettable[surveyinds]
    failedfields = targettable[failedinds]
    earlysciencefields = targettable[earlyscienceinds]
    argofields = targettable[argoinds]
    testfields = targettable[testinds]
    #get the cols want, can always/update change later if need be
    obstable = surveyfields['taskID','name','field_ra','field_dec',
                             'telescopes','duration','quality',
                             'beamPattern']
    argoobs = argofields['taskID','name','field_ra','field_dec',
                          'telescopes','duration','quality',
                          'beamPattern']
    testobs = testfields['taskID','name','field_ra','field_dec',
                          'telescopes','duration','quality',
                          'beamPattern']
    earlyobs = earlysciencefields['taskID','name','field_ra','field_dec',
                                   'telescopes','duration','quality',
                                   'beamPattern']
    failedobs = failedfields['taskID','name','field_ra','field_dec',
                             'telescopes','duration','quality',
                             'beamPattern']
    #if indicated, write tables out to a record
    #use the current date, to keep a record
    #also write to one w/ no date for regular observations (most up to date)
    #can be used by default in other things
    #get the date
    date = datetime.today().strftime('%Y-%m-%d')
    if write is True:
        ascii.write(obstable,
                    os.path.join(filedir,'obsatdb_{}.csv'.format(date)),
                    format='csv')
        ascii.write(obstable,
                    os.path.join(filedir,'obsatdb.csv'),
                    format='csv')
        ascii.write(failedobs,
                    os.path.join(filedir,'badatdb_{}.csv'.format(date)),
                    format='csv')
        ascii.write(failedobs,
                    os.path.join(filedir,'badatdb.csv'),
                    format='csv')
        ascii.write(argoobs,
                    os.path.join(filedir,'argoatdb_{}.csv'.format(date)),
                    format='csv')
        ascii.write(earlyobs,
                    os.path.join(filedir,'earlysciatdb_{}.csv'.format(date)),
                    format='csv')
        ascii.write(testobs,
                    os.path.join(filedir,'testatdb_{}.csv'.format(date)),
                    format='csv')

    #return the tables
    return obstable, argoobs, earlyobs, testobs

    

