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

def get_obslist():
    """
    Query ATDB to get observation list
    Simple but slow
    If I use this code more than once in module,
    I'm doing something wrong
    """
    obslist = atdbquery.atdbquery('imaging',False,False)
    return obslist

def get_obstable():
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
    for i,(name,taskid) in enumerate(zip(targettable['name'],targettable['taskID'])):
            #first grab specific fields that I know are tests
            if (taskid == '191030203'):
                        testinds.append(i)
                            #find ARGO fields
            elif name[0:4] == 'ARGO':
                        argoinds.append(i)
                            #find early science fields based on date
            elif (int(taskid) < 190702000) and (int(taskid) > 190409000):
                        earlyscienceinds.append(i)
                            #find survey fields based on name length and date
            elif (len(name) == 10) and (int(taskid) > 190702000):
                        surveyinds.append(i)
                            #dump everything else somewhere
            else:
                        otherinds.append(i)

    #get various table subsets
    surveyfields = targettable[surveyinds]
    otherfields = targettable[otherinds]
    earlysciencefields = targettable[earlyscienceinds]
    argofields = targettable[argoinds]
    #get the cols want, can always/update change later if need be
    obstable = surveyfields(['taskID','name','field_ra','field_dec',
                             'telescopes','duration','quality',
                             'beamPattern'])
    argoobs = argofields(['taskID','name','field_ra','field_dec',
                          'telescopes','duration','quality',
                          'beamPattern'])
    earlyobs = earlysciencefields(['taskID','name','field_ra','field_dec',
                                   'telescopes','duration','quality',
                                   'beamPattern'])
    
    #return the tables
    return obstable, argoobs, earlyobs
    #also write them out for a record
    #use the current date
    #get the date
    date = datetime.today().strftime('%Y-%m-%d')
    ascii.write(obstable,'obsatdb_{}.csv'.format(date),format='csv')
    ascii.write(argoobs,'argoatdb_{}.csv'.format(date),format='csv')
    ascii.write(earlyobs,'earlysciatdb_{}.csv'.format(date),format='csv')
    
    

