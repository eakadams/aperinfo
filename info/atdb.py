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
    #know column names from ATDB
    obstable = Table(names=['taskID','name','field_ra','field_dec',
                             'telescopes','duration','quality'])
    
    

