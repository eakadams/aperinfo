#Functions for parsing validation information

from __future__ import print_function

__author__="E.A.K. Adams"

"""
Functions related to combining
and parsing validation information.
Base files are in "files" directory
"""

import glob
import sys
from astropy.table import Table
import numpy as np
from astropy.io import ascii
from datetime import datetime
import os
import shutil
import matplotlib.pyplot as plt

#global definition (hacky) of filedir
#filedir = "../files/"
this_dir,this_filename = os.path.split(__file__)
aperinfodir = this_dir[:-4]
filedir = os.path.join(aperinfodir,"files")
#print(filedir)


#combine validation
def combine_validation():
    """
    Use knowledge of structure
    All files are local in package
    """
    #load files
    #cont & hi for now, waiting on pol
    contfile = os.path.join(filedir,"cont_allbeams.csv")
    hifile = os.path.join(filedir,"hi_allbeams.csv")
    polfile = os.path.join(filedir,"pol_allbeams.csv")
    cont = ascii.read(contfile)
    hi = ascii.read(hifile)
    pol=ascii.read(polfile)
    
    #make a combined table that has common taskids
    #find common taskids
    cont_taskids = np.unique(cont['taskid'])
    hi_taskids = np.unique(hi['Obsid'])
    pol_taskids = np.unique(pol['taskid'])
    pol_cont_taskids = np.intersect1d(cont_taskids,pol_taskids)
    common_taskids = np.intersect1d(pol_cont_taskids, hi_taskids)
    #make sure things are sorted
    common_taskids.sort()
    cont_taskids.sort()
    print("There are {} observations with continuum validation".
          format(len(cont_taskids)))
    print("There are {} observations with HI validation".
          format(len(hi_taskids)))
    print("There are {} observations with polarization validation".
          format(len(pol_taskids)))
    print("There are {} observations with all validation".
          format(len(common_taskids)))
    #setup table
    #use continuum as basis for now; other valid not up to date
    #and cont determines release
    table_length = len(common_taskids) * 40
    table_length = len(cont_taskids) * 40
    t = Table()
    t['taskid'] = np.full(table_length,'000000000')
    t['beam'] = np.empty(table_length,dtype=int)
    t['cont_pass'] = np.full(table_length,'None',dtype=object)
    t['pol_pass'] = np.full(table_length,'None',dtype=object)
    t['pol_V_pass'] = np.full(table_length,'None',dtype=object)
    t['pol_QU_pass'] = np.full(table_length,'None',dtype=object)
    t['HI_all_good'] = np.full(table_length,'None',dtype=object)
    t['HI_all_good_ok'] = np.full(table_length,'None',dtype=object)
    t['HI_c2_good_ok'] = np.full(table_length,'None',dtype=object)
    t['HI_c2_good'] = np.full(table_length,'None',dtype=object)
    t['HI_pass'] = np.full(table_length,'None',dtype=object)
    t['all_pass'] = np.full(table_length,'None',dtype=object)
    #get relevant cont metrics (used for valid)
    #so that I can plot
    t['s_in'] = np.full(table_length,np.nan)
    t['s_out'] = np.full(table_length,np.nan)
    t['rat'] = np.full(table_length,np.nan)
    t['N2'] = np.full(table_length,np.nan)
    t['Ex-2'] = np.full(table_length,np.nan)
    #get relevant pol metrics
    t['pol_s_in'] = np.full(table_length,np.nan)
    t['pol_s_out'] = np.full(table_length,np.nan)
    t['pol_rat'] = np.full(table_length,np.nan)
    t['pol_peak'] = np.full(table_length,np.nan)
    t['pol_peak_s'] = np.full(table_length,np.nan)
    t['pol_N2'] = np.full(table_length,np.nan)
    t['pol_P2'] = np.full(table_length,np.nan)
    t['pol_Ex-2'] = np.full(table_length,np.nan)
    t['pol_Ex+2'] = np.full(table_length,np.nan)
    t['pol_ftmax'] = np.full(table_length,np.nan)
    t['pol_peak_in'] = np.full(table_length,np.nan)
    t['pol_bmin'] = np.full(table_length,np.nan)
    t['pol_bmaj'] = np.full(table_length,np.nan)
    t['pol_bpa'] = np.full(table_length,np.nan)
    t['Q_bm_fg'] = np.full(table_length,np.nan)
    t['U_bm_fg'] = np.full(table_length,np.nan)
    t['Q_st_fg'] = np.full(table_length,np.nan)
    t['U_st_fg'] = np.full(table_length,np.nan)
    #get relevant HI metrics
    t['c2'] = np.full(table_length,'N')
    t['c1'] = np.full(table_length,'N')
    t['c0'] = np.full(table_length,'N')
    t['rms_c2'] = np.full(table_length,np.nan)
    t['rms_c1'] = np.full(table_length,np.nan)
    t['rms_c0'] = np.full(table_length,np.nan)
    t['lgfrac_c2'] = np.full(table_length,np.nan)
    t['lgfrac_c1'] = np.full(table_length,np.nan)
    t['lgfrac_c0'] = np.full(table_length,np.nan)
    t['prom_c2'] = np.full(table_length,np.nan)
    t['prom_c1'] = np.full(table_length,np.nan)
    t['prom_c0'] = np.full(table_length,np.nan)
    #now iterate through taskids
    #for i,tid in enumerate(common_taskids):
    for i,tid in enumerate(cont_taskids):
        #and go through the beams
        for b in range(40):
            ind = i*40 + b
            t['taskid'][ind] = tid
            t['beam'][ind] = b
            #find continuum & HI index
            cont_ind = np.where(
                (cont['taskid'] == tid) &
                (cont['beam'] == b) )[0]
            hi_ind = np.where(
                (hi['Obsid'] == tid) &
                (hi['Beam'] == b))[0]
            pol_ind =  np.where(
                (pol['taskid'] == tid) &
                (pol['beam'] == b) )[0]
            #fill in table values
            #something weird is happening with continuum so do if statement
            #t['cont_pass'][ind] = cont['pass'][cont_ind]
            #column is string, not boolean. booo.
            #need to check that index exists
            #should have length exactly one; only one index
            if len(cont_ind) == 1:
                if cont['pass'][cont_ind] == 'True':
                    t['cont_pass'][ind] = True
                else:
                    t['cont_pass'][ind] = False
                #add cont metric
                t['s_in'][ind] = cont['s_in'][cont_ind]
                t['s_out'][ind] = cont['s_out'][cont_ind]
                t['rat'][ind] = cont['rat'][cont_ind]
                t['N2'][ind] = cont['N2'][cont_ind]
                t['Ex-2'][ind] = cont['Ex-2'][cont_ind]
            #and do polarization
            if len(pol_ind) == 1:
                if pol['pass_V'][pol_ind] == 'True':
                    t['pol_V_pass'][ind] = True
                elif pol['pass_V'][pol_ind] == 'False':
                    t['pol_V_pass'][ind] = False
                
                if pol['pass_QU'][pol_ind] == 'True':
                    t['pol_QU_pass'][ind] = True
                elif pol['pass_QU'][pol_ind] == 'False':
                    t['pol_QU_pass'][ind] = False

                #do combined passing
                if ((pol['pass_V'][pol_ind] == 'True') and
                    (pol['pass_QU'][pol_ind] == 'True' )):
                    t['pol_pass'][ind] = True
                elif ((pol['pass_V'][pol_ind] == 'None') and
                      (pol['pass_QU'][pol_ind] == 'None' )):
                    t['pol_pass'][ind] = 'None'
                else:
                    t['pol_pass'][ind] = False

                #do pol metrics
                t['pol_s_in'][ind] = pol['s_in'][pol_ind]
                t['pol_s_out'][ind] = pol['s_out'][pol_ind]
                t['pol_rat'][ind] = pol['rat'][pol_ind]
                t['pol_N2'][ind] = pol['N2'][pol_ind]
                t['pol_Ex-2'][ind] = pol['Ex-2'][pol_ind]
                t['pol_ftmax'][ind] = pol['ftmax'][pol_ind]
                t['pol_peak_in'][ind] = pol['peak_in'][pol_ind]
                t['pol_peak_s'][ind] = pol['peak_s'][pol_ind]
                t['pol_peak'][ind] = pol['peak'][pol_ind]
                t['pol_P2'][ind] = pol['P2'][pol_ind]
                t['pol_bmin'][ind] = pol['bmin'][pol_ind]
                t['pol_bmaj'][ind] = pol['bmaj'][pol_ind]
                t['pol_bpa'][ind] = pol['bpa'][pol_ind]
                t['pol_Ex+2'][ind] = pol['Ex+2'][pol_ind]
                t['Q_bm_fg'][ind] = pol['Q_bm_fg'][pol_ind]
                t['U_bm_fg'][ind] = pol['U_bm_fg'][pol_ind]
                t['Q_st_fg'][ind] = pol['Q_st_fg'][pol_ind]
                t['U_st_fg'][ind] = pol['U_st_fg'][pol_ind]
    
            #for line, have to convert to boolean
            #also filling multiple columns so that I
            #can test different versions of passing
            #test all passing, that is also "pass"
            #check index exists
            if len(hi_ind) == 1:
                if hi['all_good'][hi_ind] == 1:
                    t['HI_pass'][ind] = True
                    t['HI_all_good'][ind] = True
                elif hi['all_good'][hi_ind] == 0:
                    t['HI_pass'][ind] = False
                    t['HI_all_good'][ind] = False
                #test good plus ok, all cubes
                if hi['all_good_ok'][hi_ind] == 1:
                    t['HI_all_good_ok'][ind] = True
                elif hi['all_good_ok'][hi_ind] == 0:
                    t['HI_all_good_ok'][ind] = False
                #test cube 2
                if hi['c2_good'][hi_ind] == 1:
                    t['HI_c2_good'][ind] = True
                    t['HI_c2_good_ok'][ind] = True
                elif hi['c2_ok'][hi_ind] == 1:
                    t['HI_c2_good'][ind] = False
                    t['HI_c2_good_ok'][ind] = True
                else:
                    t['HI_c2_good'][ind] = False
                    t['HI_c2_good_ok'][ind] = False

                #do hi metrics
                #cube status requires logic
                #default is bad, so just test good and okay
                if hi['c2_good'][hi_ind] == 1:
                    t['c2'][ind] = 'G'
                elif hi['c2_ok'][hi_ind] == 1:
                    t['c2'][ind] = 'O'
                else:
                    t['c2'][ind] = 'B'
                    
                if hi['c1_good'][hi_ind] == 1:
                    t['c1'][ind] = 'G'
                elif hi['c1_ok'][hi_ind] == 1:
                    t['c1'][ind] = 'O'
                else:
                    t['c1'][ind] = 'B'

                if hi['c0_good'][hi_ind] == 1:
                    t['c0'][ind] = 'G'
                elif hi['c0_ok'][hi_ind] == 1:
                    t['c0'][ind] = 'O'
                else:
                    t['c0'][ind] = 'B'

                #metrics are easy
                t['rms_c2'][ind] = hi['rms_c2'][hi_ind]
                t['rms_c1'][ind] = hi['rms_c1'][hi_ind]
                t['rms_c0'][ind] = hi['rms_c0'][hi_ind]
                t['lgfrac_c2'][ind] = hi['lgfrac_c2'][hi_ind]
                t['lgfrac_c1'][ind] = hi['lgfrac_c0'][hi_ind]
                t['lgfrac_c0'][ind] = hi['lgfrac_c1'][hi_ind]
                t['prom_c2'][ind] = hi['prom_c2'][hi_ind]
                t['prom_c1'][ind] = hi['prom_c0'][hi_ind]
                t['prom_c0'][ind] = hi['prom_c1'][hi_ind]



    #now print and return some useful information
    pass_cont = np.where(t['cont_pass'] == True)[0]
    pass_hi = np.where(t['HI_pass'] == True)[0]
    pass_both = np.where(
        (t['cont_pass'] == True) &
        (t['HI_pass'] == True))[0]
    print("There are {} beams in total".format(len(t)))
    print("{} beams pass continuum validation".format(len(pass_cont)))
    print("{} beams pass HI validation".format(len(pass_hi)))
    print("{} beams pass both HI and cont validation".format(len(pass_both)))

    #add in polarization
    pass_V_pol = np.where(t['pol_V_pass'] == True)[0]
    pass_QU_pol = np.where(t['pol_QU_pass'] == True)[0]

    pass_all = np.where(
        (t['cont_pass'] == True) &
        (t['HI_pass'] == True) &
        (t['pol_V_pass'] == True) &
        (t['pol_QU_pass'] == True) )[0]

    print("{} beams passed all  validation".format(len(pass_all)))

    #save table to file
    ascii.write(t,
                os.path.join(filedir,'combined_valid.csv'),
                overwrite = True,format='csv')

def compare_hi_cont_valid():
    """
    Compare HI and continuum validation
    """
    #read file in
    t = ascii.read(os.path.join(filedir,'combined_valid.csv'))
    pass_cont = np.where(t['cont_pass'] == 'True')[0]
    pass_hi = np.where(t['HI_pass'] == 'True')[0]
    pass_both = np.where(
        (t['cont_pass'] == 'True') &
        (t['HI_pass'] == 'True'))[0]
    #want to get those that pass continuum and not line
    #so they can be further investigated
    #and maybe also the other way around
    #and then maybe it's interesting to plot those points also?
    pass_cont_not_hi = np.where(
        (t['cont_pass'] == 'True') &
        (t['HI_pass'] == 'False'))[0]
    pass_hi_not_cont =  np.where(
        (t['cont_pass'] == 'False') &
        (t['HI_pass'] == 'True'))[0]

    #do some testing for other types of HI passing
    pass_hi_ok = np.where(t['HI_all_good_ok'] == 'True')[0]
    pass_both_ok = np.where(
        (t['cont_pass'] == 'True') &
        (t['HI_all_good_ok'] == 'True'))[0]
    pass_hi_c2_ok = np.where(t['HI_c2_good_ok'] == 'True')[0]
    pass_both_c2_ok = np.where(
        (t['cont_pass'] == 'True') &
        (t['HI_c2_good_ok'] == 'True'))[0]
    pass_hi_c2 = np.where(t['HI_c2_good'] == 'True')[0]
    pass_both_c2 = np.where(
        (t['cont_pass'] == 'True') &
        (t['HI_c2_good'] == 'True'))[0]

    print("{} beams pass HI good/ok".format(len(pass_hi_ok)))
    print("{} beams pass both HI good/ok and cont".format(len(pass_both_ok)))

    print(("{} more beams passed HI ok than just good. "
           "{} more beams passed cont plus HI then").
          format((len(pass_hi_ok)-len(pass_hi)),(len(pass_both_ok)-len(pass_both))))

    print("{} beams pass HI c2".format(len(pass_hi_c2)))
    print("{} beams pass both HI c2 and cont".format(len(pass_both_c2)))
    
    print("{} beams pass HI c2 good/ok".format(len(pass_hi_c2_ok)))
    print("{} beams pass both HI c2 good/ok and cont".format(len(pass_both_c2_ok)))


    #make some plots of continuum metrics
    #first inner vs. outer noise
    fig, ax = plt.subplots(figsize=(8,8))
    ax.scatter(t['s_in'],t['s_out'],marker='.',s=5,color='grey',label='All beams')
    ax.scatter(t['s_in'][pass_cont],t['s_out'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    ax.scatter(t['s_in'][pass_hi],t['s_out'][pass_hi],marker='s',
               s=20,label='Good HI',
               facecolors='none',edgecolors='#ff7f0e')
    ax.set_xlabel('Inner noise')
    ax.set_ylabel('Outer noise')
    ax.set_xlim([25,100])
    ax.set_ylim([25,100])
    ax.plot([60,60],[0,1000],'k:')
    ax.plot([0,1000],[60,60],'k:')
    plt.legend()
    pltname = 'inner_outer_noise.png'
    plt.savefig(os.path.join(filedir,pltname))

    #histogram of inner and outer noise
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12,6))
    binvals = np.arange(25,100,0.1)
    ax1.hist([t['s_in'],t['s_in'][pass_cont],t['s_in'][pass_hi]],
             bins=binvals,
             label=['All','Valid cont', 'Good HI'])
    ax1.set_xlabel('Inner noise')
    ax2.hist([t['s_out'],t['s_out'][pass_cont],t['s_out'][pass_hi]],
             bins=binvals,
             label=['All','Valid cont', 'Good HI'])
    ax2.set_xlabel('Outer noise')
    plt.legend()
    pltname = 'noise_histogram.png'
    plt.savefig(os.path.join(filedir,pltname))

    #histogram of other criteria
    fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(18,6))
    ratbins = np.arange(1,4,0.1)
    ax1.hist([t['rat'],t['rat'][pass_cont],t['rat'][pass_hi]],
             bins=ratbins,
             label=['All','Valid cont', 'Good HI'])
    ax1.set_xlabel('Ratio inner-to-outer noise')
    n2bins = np.arange(1,15,.1)
    ax2.hist([t['N2'],t['N2'][pass_cont],t['N2'][pass_hi]],
             bins=n2bins,
             label=['All','Valid cont', 'Good HI'])
    ax2.set_xlabel('N2')
    ex2bins = np.arange(0,4000,10)
    ax3.hist([t['Ex-2'],t['Ex-2'][pass_cont],t['Ex-2'][pass_hi]],
             bins=ex2bins,
             label=['All','Valid cont', 'Good HI'])
    ax3.set_xlabel('Ex-2')
    ax1.legend()
    pltname = 'hist_other_cont.png'
    plt.savefig(os.path.join(filedir,pltname))

    #compare other metrics
    fig, ((ax1, ax2, ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(18,12))
    #ratio vs N2
    ax1.scatter(t['rat'],t['N2'],marker='.',s=5,color='grey',label='All beams')
    ax1.scatter(t['rat'][pass_cont],t['N2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax1.scatter(t['rat'][pass_hi],t['N2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax1.scatter(t['rat'][pass_hi_not_cont],t['N2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax1.scatter(t['rat'][pass_cont_not_hi],t['N2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')
    
    ax1.set_xlabel('Ratio inner-to-outer noise')
    ax1.set_ylabel('N2')
    ax1.set_xlim([0.8,4])
    ax1.set_ylim([0.8,10])
    ax1.plot([1.225,1.225],[0,1000],'k--')
    ax1.plot([1.15,1.15],[0,1000],'k:')
    ax1.plot([0,1000],[4.5,4.5],'k:')
    #zoom in on valid area
    ax4.scatter(t['rat'],t['N2'],marker='.',s=5,color='grey',label='All beams')
    ax4.scatter(t['rat'][pass_cont],t['N2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax4.scatter(t['rat'][pass_hi],t['N2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax4.scatter(t['rat'][pass_hi_not_cont],t['N2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax4.scatter(t['rat'][pass_cont_not_hi],t['N2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')

    ax4.set_xlabel('Ratio inner-to-outer noise')
    ax4.set_ylabel('N2')
    ax4.set_xlim([1.0,2.0])
    ax4.set_ylim([0.8,4.8])
    ax4.plot([1.225,1.225],[0,1000],'k--')
    ax4.plot([1.15,1.15],[0,1000],'k:')
    ax4.plot([0,1000],[4.5,4.5],'k:')

    #ratio vs Ex-2
    print(np.nanmax(t['Ex-2']))
    ax2.scatter(t['rat'],t['Ex-2'],marker='.',s=5,color='grey',label='All beams')
    ax2.scatter(t['rat'][pass_cont],t['Ex-2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax2.scatter(t['rat'][pass_hi],t['Ex-2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax2.scatter(t['rat'][pass_hi_not_cont],t['Ex-2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax2.scatter(t['rat'][pass_cont_not_hi],t['Ex-2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')
    ax2.set_xlabel('Ratio inner-to-outer noise')
    ax2.set_ylabel('Ex-2')
    ax2.set_xlim([0.8,4])
    ax2.set_ylim([0,5050])
    ax2.plot([1.225,1.225],[0,5050],'k--')
    ax2.plot([1.15,1.15],[0,5050],'k:')
    ax2.plot([0,1000],[400,400],'k:')
    #zoom in on valid area
    ax5.scatter(t['rat'],t['Ex-2'],marker='.',s=5,color='grey',label='All beams')
    ax5.scatter(t['rat'][pass_cont],t['Ex-2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax5.scatter(t['rat'][pass_hi],t['Ex-2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax5.scatter(t['rat'][pass_hi_not_cont],t['Ex-2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax5.scatter(t['rat'][pass_cont_not_hi],t['Ex-2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')
    ax5.set_xlabel('Ratio inner-to-outer noise')
    ax5.set_ylabel('Ex-2')
    ax5.set_xlim([1.0,2.0])
    ax5.set_ylim([0,1000])
    ax5.plot([1.225,1.225],[0,5050],'k--')
    ax5.plot([1.15,1.15],[0,5050],'k:')
    ax5.plot([0,1000],[400,400],'k:')

    #N2 vs Ex-2
    ax3.scatter(t['N2'],t['Ex-2'],marker='.',s=5,color='grey',label='All beams')
    ax3.scatter(t['N2'][pass_cont],t['Ex-2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax3.scatter(t['N2'][pass_hi],t['Ex-2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax3.scatter(t['N2'][pass_hi_not_cont],t['Ex-2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax3.scatter(t['N2'][pass_cont_not_hi],t['Ex-2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')
    ax3.set_xlabel('N2')
    ax3.set_ylabel('Ex-2')
    ax3.set_xlim([0.8,10])
    ax3.set_ylim([0,5050])
    ax3.plot([4.5,4.5],[0,10000],'k:')
    ax3.plot([0,1000],[400,400],'k:')
    #zoom in on valid area
    ax6.scatter(t['N2'],t['Ex-2'],marker='.',s=5,color='grey',label='All beams')
    ax6.scatter(t['N2'][pass_cont],t['Ex-2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax6.scatter(t['N2'][pass_hi],t['Ex-2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax6.scatter(t['N2'][pass_hi_not_cont],t['Ex-2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax6.scatter(t['N2'][pass_cont_not_hi],t['Ex-2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')
    ax6.set_xlabel('N2')
    ax6.set_ylabel('Ex-2')
    ax6.set_xlim([0.8,4.8])
    ax6.set_ylim([0,1000])
    ax6.plot([4.5,4.5],[0,5050],'k:')
    ax6.plot([0,1000],[400,400],'k:')

    ax1.legend()
    pltname = 'other_cont_crit.png'
    plt.savefig(os.path.join(filedir,pltname))

    ascii.write(t[pass_cont_not_hi],
                os.path.join(filedir,'combined_valid_pass_cont_not_hi.csv'),
                overwrite = True,format='csv')

    
def compare_pol_cont_valid():
    """
    Compare continuum and polarization validation
    """
    #read file in
    t = ascii.read(os.path.join(filedir,'combined_valid.csv'))
    pass_cont = np.where(t['cont_pass'] == 'True')[0]
    pass_pol_V = np.where(t['pol_V_pass'] == 'True')[0]
    pass_pol_QU = np.where(t['pol_QU_pass'] == 'True')[0]
    pass_both = np.where(
        (t['cont_pass'] == 'True') &
        (t['pol_V_pass'] == 'True') &
        (t['pol_QU_pass'] == 'True') )[0]
    #want to get those that pass continuum and not pol
    #so they can be further investigated
    #and maybe also the other way around
    #and then maybe it's interesting to plot those points also?
    #also have to think about which pol is passed
    #actually, first just focus on QU vs V validation
    pass_V_not_QU = np.where(
        (t['pol_V_pass'] == 'True') &
        (t['pol_QU_pass'] == 'False') )[0]
    pass_QU_not_V = np.where(
        (t['pol_QU_pass'] == 'True') &
        (t['pol_V_pass'] == 'False') )[0]
    
    pass_cont_not_V = np.where(
        (t['cont_pass'] == 'True') &
        (t['pol_V_pass'] == 'False'))[0]
    pass_cont_not_QU =  np.where(
        (t['cont_pass'] == 'True') &
        (t['pol_QU_pass'] == 'False'))[0]

    

    print("{} beams pass Stokes V".format(len(pass_pol_V)))
    print("{} beams pass Stokes QU".format(len(pass_pol_QU)))
    print("{} beams pass continuum".format(len(pass_cont)))
    print("{} beams pass continuum and (all) Stokes".format(len(pass_both)))
    print("{} beams pass continuum and fail (all) pol".format(len(pass_cont)-len(pass_both)))


    #first how many pass V but not QU?
    #that is not quite answered by printout above

    pass_V_not_QU = np.where(
        (t['pol_V_pass'] == 'True') &
        (t['pol_QU_pass'] == 'False'))[0]
    #print("{} beams pass Stokes V but not QU".format(len(pass_V_not_QU)))

    #then the question is whether these pass continuum or not
    #that tells me whether I have to worry about this or not
    #want to find pass cont, V but not QU
    array_pass_cont_and_V = np.empty(len(t),dtype=bool)
    for i in enumerate(array_pass_cont_and_V):
        if (t['pol_V_pass'][i] == 'True') and (t['cont_pass'][i] == 'True'):
            array_pass_cont_and_V[i] = True
        else:
            array_pass_cont_and_V[i] = False
    pass_contV_not_QU = np.where( (array_pass_cont_and_V == True) &
                                  (t['pol_QU_pass'] == 'False') )[0]
    
    print("{} beams pass continuum and Stokes V but not QU".format(len(pass_contV_not_QU)))
    #so there are 209 beams to potentially worry about here
    #this is about ~50% of the beams that pass continuum but not pol

    #first check beam issue - this shouldn't be it
    #check for where pass cont and Q_bm_fg > .3
    pass_cont_notQbeam = np.where( (t['cont_pass'] == 'True') &
                                   (t['Q_bm_fg'] > 0.33) )[0]
    pass_contV_notQbeam = np.where( (array_pass_cont_and_V == True) &
                                   (t['Q_bm_fg'] > 0.33) )[0]

    print("{} beams pass continuum but not Q beam".format(len(pass_cont_notQbeam)))
    print("{} beams pass continuum/V but not Q beam".format(len(pass_contV_notQbeam)))

    #print these out because they're possibly interesting
    ascii.write(t[pass_cont_notQbeam],
                os.path.join(filedir,'pass_cont_notQbeam.csv'),
                overwrite = True,format='csv')

    #so, about half that fail polarization are becasue of QU beam size
    #now investigate those which pass continuum but fail V

    print("{} beams pass continuum but not V".format(len(pass_cont_not_V)))
    #find which V metrics they fail
    fail_sigma_in = np.where(t['pol_s_in'][pass_cont_not_V] > 60.)[0]
    fail_sigma_out = np.where(t['pol_s_out'][pass_cont_not_V] > 60.)[0]

    print("{} beams fail because of sigma in".format(len(fail_sigma_in)))
    print("{} beams fail because of sigma out".format(len(fail_sigma_out)))
    #not a major cause

    fail_ftmax = np.where(t['pol_ftmax'][pass_cont_not_V] > 25)[0]
    fail_pin = np.where(t['pol_peak_in'][pass_cont_not_V] > 4000)[0]
    print("{} beams fail because of FTmax".format(len(fail_ftmax)))
    print("{} beams fail because of peak inner".format(len(fail_pin)))
    #peak inner seems to be main cause - all those that fail V are because of it
    #doesn't mean other things might not contribute

    print(np.nanmin(t['pol_peak_in'][pass_cont_not_V]),np.nanmax(t['pol_peak_in'][pass_cont_not_V]),np.nanmedian(t['pol_peak_in'][pass_cont_not_V]))

    #but what are units?
    #look at this more closely, make histograms.....
    #think that if I treat peak_in as uJy, then failures are no longer so bad.....

    fail_rat = np.where( (t['pol_rat'][pass_cont_not_V] < 0.99) &
                         (t['pol_rat'][pass_cont_not_V] > 1.1 ) )[0]
    fail_combi = np.where( (t['pol_rat'][pass_cont_not_V] > 1.05) &
                           (t['pol_N2'][pass_cont_not_V] > 1.1) &
                           (t['pol_Ex-2'][pass_cont_not_V] > 1000) &
                           (t['pol_P2'][pass_cont_not_V] > 1.1 ) )[0]
    fail_beam = np.where( t['pol_bmin'][pass_cont_not_V] > 15. )[0]
    fail_n2_p2 = np.where( (t['pol_N2'][pass_cont_not_V] < 0.9) &
                           (t['pol_P2'][pass_cont_not_V] < 0.9) )[0]
    
    print("{} beams fail because of ratio requirements".format(len(fail_rat)))
    print("{} beams fail because of combination requirement".format(len(fail_combi)))
    print("{} beams fail because of beam requirements".format(len(fail_beam)))
    print("{} beams fail because of N2/P2 requirements".format(len(fail_n2_p2)))

    #combined requirement is one that causest most beams to fail
    #if changed tomatch tom, what happens?
    fail_combi_tom = np.where( (t['pol_rat'][pass_cont_not_V] > 1.15) &
                               (t['pol_N2'][pass_cont_not_V] > 4.5) &
                               (t['pol_Ex-2'][pass_cont_not_V] > 400) )[0]
    
    print("{} beams fail when same combination requirement as cont is used".format(len(fail_combi_tom)))

    
    """ 
    #make some plots of continuum metrics
    #first inner vs. outer noise
    fig, ax = plt.subplots(figsize=(8,8))
    ax.scatter(t['s_in'],t['s_out'],marker='.',s=5,color='grey',label='All beams')
    ax.scatter(t['s_in'][pass_cont],t['s_out'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    ax.scatter(t['s_in'][pass_hi],t['s_out'][pass_hi],marker='s',
               s=20,label='Good HI',
               facecolors='none',edgecolors='#ff7f0e')
    ax.set_xlabel('Inner noise')
    ax.set_ylabel('Outer noise')
    ax.set_xlim([25,100])
    ax.set_ylim([25,100])
    ax.plot([60,60],[0,1000],'k:')
    ax.plot([0,1000],[60,60],'k:')
    plt.legend()
    pltname = 'inner_outer_noise.png'
    plt.savefig(os.path.join(filedir,pltname))

    #histogram of inner and outer noise
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12,6))
    binvals = np.arange(25,100,0.1)
    ax1.hist([t['s_in'],t['s_in'][pass_cont],t['s_in'][pass_hi]],
             bins=binvals,
             label=['All','Valid cont', 'Good HI'])
    ax1.set_xlabel('Inner noise')
    ax2.hist([t['s_out'],t['s_out'][pass_cont],t['s_out'][pass_hi]],
             bins=binvals,
             label=['All','Valid cont', 'Good HI'])
    ax2.set_xlabel('Outer noise')
    plt.legend()
    pltname = 'noise_histogram.png'
    plt.savefig(os.path.join(filedir,pltname))

    #histogram of other criteria
    fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(18,6))
    ratbins = np.arange(1,4,0.1)
    ax1.hist([t['rat'],t['rat'][pass_cont],t['rat'][pass_hi]],
             bins=ratbins,
             label=['All','Valid cont', 'Good HI'])
    ax1.set_xlabel('Ratio inner-to-outer noise')
    n2bins = np.arange(1,15,.1)
    ax2.hist([t['N2'],t['N2'][pass_cont],t['N2'][pass_hi]],
             bins=n2bins,
             label=['All','Valid cont', 'Good HI'])
    ax2.set_xlabel('N2')
    ex2bins = np.arange(0,4000,10)
    ax3.hist([t['Ex-2'],t['Ex-2'][pass_cont],t['Ex-2'][pass_hi]],
             bins=ex2bins,
             label=['All','Valid cont', 'Good HI'])
    ax3.set_xlabel('Ex-2')
    ax1.legend()
    pltname = 'hist_other_cont.png'
    plt.savefig(os.path.join(filedir,pltname))

    #compare other metrics
    fig, ((ax1, ax2, ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(18,12))
    #ratio vs N2
    ax1.scatter(t['rat'],t['N2'],marker='.',s=5,color='grey',label='All beams')
    ax1.scatter(t['rat'][pass_cont],t['N2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax1.scatter(t['rat'][pass_hi],t['N2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax1.scatter(t['rat'][pass_hi_not_cont],t['N2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax1.scatter(t['rat'][pass_cont_not_hi],t['N2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')
    
    ax1.set_xlabel('Ratio inner-to-outer noise')
    ax1.set_ylabel('N2')
    ax1.set_xlim([0.8,4])
    ax1.set_ylim([0.8,10])
    ax1.plot([1.225,1.225],[0,1000],'k--')
    ax1.plot([1.15,1.15],[0,1000],'k:')
    ax1.plot([0,1000],[4.5,4.5],'k:')
    #zoom in on valid area
    ax4.scatter(t['rat'],t['N2'],marker='.',s=5,color='grey',label='All beams')
    ax4.scatter(t['rat'][pass_cont],t['N2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax4.scatter(t['rat'][pass_hi],t['N2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax4.scatter(t['rat'][pass_hi_not_cont],t['N2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax4.scatter(t['rat'][pass_cont_not_hi],t['N2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')

    ax4.set_xlabel('Ratio inner-to-outer noise')
    ax4.set_ylabel('N2')
    ax4.set_xlim([1.0,2.0])
    ax4.set_ylim([0.8,4.8])
    ax4.plot([1.225,1.225],[0,1000],'k--')
    ax4.plot([1.15,1.15],[0,1000],'k:')
    ax4.plot([0,1000],[4.5,4.5],'k:')

    #ratio vs Ex-2
    print(np.nanmax(t['Ex-2']))
    ax2.scatter(t['rat'],t['Ex-2'],marker='.',s=5,color='grey',label='All beams')
    ax2.scatter(t['rat'][pass_cont],t['Ex-2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax2.scatter(t['rat'][pass_hi],t['Ex-2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax2.scatter(t['rat'][pass_hi_not_cont],t['Ex-2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax2.scatter(t['rat'][pass_cont_not_hi],t['Ex-2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')
    ax2.set_xlabel('Ratio inner-to-outer noise')
    ax2.set_ylabel('Ex-2')
    ax2.set_xlim([0.8,4])
    ax2.set_ylim([0,5050])
    ax2.plot([1.225,1.225],[0,5050],'k--')
    ax2.plot([1.15,1.15],[0,5050],'k:')
    ax2.plot([0,1000],[400,400],'k:')
    #zoom in on valid area
    ax5.scatter(t['rat'],t['Ex-2'],marker='.',s=5,color='grey',label='All beams')
    ax5.scatter(t['rat'][pass_cont],t['Ex-2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax5.scatter(t['rat'][pass_hi],t['Ex-2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax5.scatter(t['rat'][pass_hi_not_cont],t['Ex-2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax5.scatter(t['rat'][pass_cont_not_hi],t['Ex-2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')
    ax5.set_xlabel('Ratio inner-to-outer noise')
    ax5.set_ylabel('Ex-2')
    ax5.set_xlim([1.0,2.0])
    ax5.set_ylim([0,1000])
    ax5.plot([1.225,1.225],[0,5050],'k--')
    ax5.plot([1.15,1.15],[0,5050],'k:')
    ax5.plot([0,1000],[400,400],'k:')

    #N2 vs Ex-2
    ax3.scatter(t['N2'],t['Ex-2'],marker='.',s=5,color='grey',label='All beams')
    ax3.scatter(t['N2'][pass_cont],t['Ex-2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax3.scatter(t['N2'][pass_hi],t['Ex-2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax3.scatter(t['N2'][pass_hi_not_cont],t['Ex-2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax3.scatter(t['N2'][pass_cont_not_hi],t['Ex-2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')
    ax3.set_xlabel('N2')
    ax3.set_ylabel('Ex-2')
    ax3.set_xlim([0.8,10])
    ax3.set_ylim([0,5050])
    ax3.plot([4.5,4.5],[0,10000],'k:')
    ax3.plot([0,1000],[400,400],'k:')
    #zoom in on valid area
    ax6.scatter(t['N2'],t['Ex-2'],marker='.',s=5,color='grey',label='All beams')
    ax6.scatter(t['N2'][pass_cont],t['Ex-2'][pass_cont],marker='o',
               s=10,label='Valid continuum')
    #ax6.scatter(t['N2'][pass_hi],t['Ex-2'][pass_hi],marker='s',
    #           s=20,label='Good HI',
    #           facecolors='none',edgecolors='#ff7f0e')
    ax6.scatter(t['N2'][pass_hi_not_cont],t['Ex-2'][pass_hi_not_cont],marker='o',
               s=10,label='Valid HI, fail cont')
    ax6.scatter(t['N2'][pass_cont_not_hi],t['Ex-2'][pass_cont_not_hi],marker='o',
               s=10,label='Valid continuum, fail HI')
    ax6.set_xlabel('N2')
    ax6.set_ylabel('Ex-2')
    ax6.set_xlim([0.8,4.8])
    ax6.set_ylim([0,1000])
    ax6.plot([4.5,4.5],[0,5050],'k:')
    ax6.plot([0,1000],[400,400],'k:')

    ax1.legend()
    pltname = 'other_cont_crit.png'
    plt.savefig(os.path.join(filedir,pltname))

   

    ascii.write(t[pass_cont_not_hi],
                os.path.join(filedir,'combined_valid_pass_cont_not_hi.csv'),
                overwrite = True,format='csv')
    """
    
    
#combine continuum informaiton
def combine_continuum():
    """
    Use knowledge of structure
    """
    #cont valid directory
    contdir = os.path.join(filedir,"cont_valid")
    #get taskid list; this is only taskids, not full paths
    #this is a list of paths
    taskdirlist = glob.glob(
        os.path.join(contdir,
                     "[1-2][0-9][0-1][0-9][0-3][0-9][0-9][0-9][0-9]"))
    #taskidlist = os.listdir(contdir)
    taskdirlist.sort()
    #setup table that I need
    #40 * n_tid for length
    table_length = len(taskdirlist)*40
    t = Table()
    t['taskid'] = np.full(table_length,'000000000')
    t['beam'] = np.empty(table_length,dtype=int)
    t['pass'] = np.full(table_length,'None',dtype=object)
    t['s_in'] = np.full(table_length,np.nan)
    t['s_out'] = np.full(table_length,np.nan) 
    t['rat'] = np.full(table_length,np.nan) 
    t['peak_mJy'] = np.full(table_length,np.nan) 
    t['peak_sig'] = np.full(table_length,np.nan) 
    t['rusc-'] = np.full(table_length,np.nan) 
    t['rusc+'] = np.full(table_length,np.nan) 
    t['D'] = np.full(table_length,np.nan) 
    t['N2'] = np.full(table_length,np.nan) 
    t['P2'] = np.full(table_length,np.nan) 
    t['Ex-2'] = np.full(table_length,np.nan) 
    t['Ex+2'] = np.full(table_length,np.nan) 

    #now iterate through taskids
    for i,taskdir in enumerate(taskdirlist):
        #read file in
        valid_file = os.path.join(taskdir,"dynamicRange.dat")
        valid = ascii.read(valid_file)
        #now i need to fill everything
        #but have to worry about missing beams
        #taskid covers full range i*40,(i+1)*40 -1
        for b in range(40):
            ind = i*40 + b
            taskid = taskdir[-9:]
            t['taskid'][ind] = taskid
            t['beam'][ind] = b
            valid_ind = np.where(valid['col1'] == b)[0]
            #if there's a matching beam, fill information
            if len(valid_ind) == 1:
                #check for pass as True/Fail
                if valid['col2'][valid_ind] == '.':
                    t['pass'][ind] = True
                elif (valid['col2'][valid_ind] == 'X'):
                    t['pass'][ind] = False
                #fill in rest of columns
                t['s_in'][ind] = valid['col3'][valid_ind]
                t['s_out'][ind] = valid['col4'][valid_ind]
                t['rat'][ind] = valid['col5'][valid_ind]
                t['peak_mJy'][ind] = valid['col6'][valid_ind]
                t['peak_sig'][ind] = valid['col7'][valid_ind]
                t['rusc-'][ind] = valid['col8'][valid_ind]
                t['rusc+'][ind] = valid['col9'][valid_ind]
                t['D'][ind] = valid['col10'][valid_ind]
                t['N2'][ind] = valid['col11'][valid_ind]
                t['P2'][ind] = valid['col12'][valid_ind]
                t['Ex-2'][ind] = valid['col13'][valid_ind]
                t['Ex+2'][ind] = valid['col14'][valid_ind]
            else:
                #make sure pass is false if beam doesn't exist
                t['pass'][ind] = False

    #at end, write table out to csv
    ascii.write(t,
                os.path.join(filedir,'cont_allbeams.csv'),
                overwrite = True,format='csv')



#combine polarization informaiton
def combine_pol():
    """
    Use knowledge of structure
    Note that polarization informaiton is updated via:
    rsync -arvt happili-05.astron.nl:/tank/apertif/qa/polarisation/ pol_valid/
    """
    #pol valid directory
    poldir = os.path.join(filedir,"pol_valid")
    #get taskid list; this is only taskids, not full paths
    #this is a list of paths
    taskdirlist = glob.glob(
        os.path.join(poldir,
                     "[1-2][0-9][0-1][0-9][0-3][0-9][0-9][0-9][0-9]"))
    #taskidlist = os.listdir(contdir)
    taskdirlist.sort()
    #setup table that I need
    #40 * n_tid for length
    table_length = len(taskdirlist)*40
    t = Table()
    t['taskid'] = np.full(table_length,'000000000')
    t['beam'] = np.empty(table_length,dtype=int)
    t['pass_V'] = np.full(table_length,'None',dtype=object)
    t['pass_QU'] = np.full(table_length,'None',dtype=object)
    t['Q_bm_fg'] = np.full(table_length,np.nan)
    t['U_bm_fg'] = np.full(table_length,np.nan) 
    t['Q_st_fg'] = np.full(table_length,np.nan) 
    t['U_st_fg'] = np.full(table_length,np.nan) 
    t['s_in'] = np.full(table_length,np.nan) 
    t['s_out'] = np.full(table_length,np.nan) 
    t['rat'] = np.full(table_length,np.nan) 
    t['peak'] = np.full(table_length,np.nan)
    t['peak_s'] = np.full(table_length,np.nan)
    t['peak_in'] = np.full(table_length,np.nan) 
    t['N2'] = np.full(table_length,np.nan) 
    t['P2'] = np.full(table_length,np.nan) 
    t['Ex-2'] = np.full(table_length,np.nan) 
    t['Ex+2'] = np.full(table_length,np.nan)
    t['ftmax'] = np.full(table_length,np.nan)
    t['bmaj'] = np.full(table_length,np.nan)
    t['bmin'] = np.full(table_length,np.nan)
    t['bpa'] = np.full(table_length,np.nan) 

    #empty table as plaveholder for missing files
    foo = Table()
    foo['B'] = np.full(1,np.nan)
    
    #now iterate through taskids
    #have two files for every one - bleh
    for i,taskdir in enumerate(taskdirlist):
        #read file in
        #have to add a bunch of try/excepts because there are missing files
        #why?!?!
        valid_file_V = os.path.join(taskdir,"dynamicRange_V.dat")
        try:
            valid_V = ascii.read(valid_file_V,header_start=4)
        except:
            print("{} doesn't exist".format(valid_file_V))
            valid_V = foo
        valid_file_QU = os.path.join(taskdir,"dynamicRange_QU.dat")
        try:
            valid_QU = ascii.read(valid_file_QU,header_start=4)
        except:
            print("{} doesn't exist".format(valid_file_QU))
            valid_QU = foo
        #now i need to fill everything
        #but have to worry about missing beams
        #taskid covers full range i*40,(i+1)*40 -1
        for b in range(40):
            ind = i*40 + b
            taskid = taskdir[-9:]
            t['taskid'][ind] = taskid
            t['beam'][ind] = b
            #think in this case there is always information
            #but still check for each valid
            valid_ind_V = np.where(valid_V['B'] == b)[0]
            valid_ind_QU = np.where(valid_QU['B'] == b)[0]
            #if there's a matching beam, fill information
            if len(valid_ind_V) == 1:
                #check for pass as True/Fail
                if valid_V['Q'][valid_ind_V] == '.':
                    t['pass_V'][ind] = True
                elif ( (valid_V['Q'][valid_ind_V] == 'X') and
                       np.isnan(valid_V['s_in'][valid_ind_V]) ):
                    t['pass_V'][ind] = 'None'
                    #print(taskid,b)
                elif valid_V['Q'][valid_ind_V] == 'X':
                    t['pass_V'][ind] = False
                #fill in rest of columns
                t['s_in'][ind] = valid_V['s_in'][valid_ind_V]
                t['s_out'][ind] = valid_V['s_out'][valid_ind_V]
                t['rat'][ind] = valid_V['rat'][valid_ind_V]
                t['peak'][ind] = valid_V['peak'][valid_ind_V]
                t['peak_s'][ind] = valid_V['peak_s'][valid_ind_V]
                t['peak_in'][ind] = valid_V['peak_in'][valid_ind_V]
                t['N2'][ind] = valid_V['N2'][valid_ind_V]
                t['P2'][ind] = valid_V['P2'][valid_ind_V]
                t['Ex-2'][ind] = valid_V['Ex-2'][valid_ind_V]
                t['Ex+2'][ind] = valid_V['Ex+2'][valid_ind_V]
                t['ftmax'][ind] = valid_V['ftmax'][valid_ind_V]
                t['bmaj'][ind] = valid_V['bmaj'][valid_ind_V]
                t['bmin'][ind] = valid_V['bmin'][valid_ind_V]
                t['bpa'][ind] = valid_V['bpa'][valid_ind_V]
            else:
                #make sure pass is false if beam doesn't exist
                t['pass_V'][ind] = False
            #repeat for QU
            if len(valid_ind_QU) == 1:
                #check for pass as True/Fail
                if valid_QU['Q'][valid_ind_QU] == '.':
                    t['pass_QU'][ind] = True
                elif ( (valid_QU['Q'][valid_ind_QU] == 'X') and
                       (valid_QU['Q_bm_fg'][valid_ind_QU] == 1.0) ):
                    t['pass_QU'][ind] = 'None'
                elif valid_QU['Q'][valid_ind_QU] == 'X':
                    t['pass_QU'][ind] = False
                #fill in rest of columns
                t['Q_bm_fg'][ind] = valid_QU['Q_bm_fg'][valid_ind_QU]
                t['U_bm_fg'][ind] = valid_QU['U_bm_fg'][valid_ind_QU]
                t['Q_st_fg'][ind] = valid_QU['Q_st_fg'][valid_ind_QU]
                t['U_st_fg'][ind] = valid_QU['U_st_fg'][valid_ind_QU]
            else:
                #make sure pass is false if beam doesn't exist
                t['pass_QU'][ind] = False
    #at end, write table out to csv
    ascii.write(t,
                os.path.join(filedir,'pol_allbeams.csv'),
                overwrite = True,format='csv')

def do_pol_valid():
    """
    Use the local csv file to update the polarization validation
    This can account for changing criteria
    Motivated 30 Sep 2020 by change from 12.5" to 15" for b_min
    """
    poltable = ascii.read(os.path.join(filedir,'pol_allbeams.csv'))
    QU_pass_new = np.full(len(poltable),True,dtype=object)
    V_pass_new = np.full(len(poltable),True,dtype=object)
    #find where V fails
    failind_sin = np.where(poltable['s_in'] > 60.)[0]
    failind_sout = np.where(poltable['s_out'] > 60.)[0]
    failind_bmin = np.where(poltable['bmin'] > 15.)[0]
    failind_ftmax = np.where(poltable['ftmax'] > 25.)[0]
    failind_peak_in = np.where(poltable['peak_in'] > 4000.)[0]
    failind_nodata = np.hstack(np.argwhere(np.isnan(poltable['s_in'])))
    failind_V = np.unique(np.concatenate((failind_sin,failind_sout,failind_bmin,
                               failind_ftmax,failind_peak_in,failind_nodata)))

    V_pass_new[failind_V] = False
    noind_V = np.where(poltable['pass_V'] == 'None')[0]
    V_pass_new[noind_V] = 'None'
    poltable['pass_V'] = V_pass_new
    
    #check QU also, but should be good

    #now check for where things fails
    failind_Qbm = np.where(poltable['Q_bm_fg'] > 0.31)[0]
    failind_Ubm = np.where(poltable['U_bm_fg'] > 0.31)[0]
    failind_Qst = np.where(poltable['Q_st_fg'] > 0.31)[0]
    failind_Ust = np.where(poltable['U_st_fg'] > 0.31)[0]

    failind_QU = np.unique(np.concatenate((failind_Qbm,failind_Ubm,
                                           failind_Qst,failind_Ust)))

    #nothing actually changed, as I hoped/expected
    QU_pass_new[failind_QU] = False
    noind_QU = np.where(poltable['pass_QU'] == 'None')[0]
    QU_pass_new[noind_QU] = 'None'
    poltable['pass_QU'] = QU_pass_new

    #update to set 200309042 as nan
    #want to do this quickly, in bulk, but doesn't seem possible
    badtask = np.where(poltable['taskid'] == 200309042)[0]
    print(poltable[badtask])
    cols = poltable.colnames
    poltable[cols[2]][badtask] = np.full(len(badtask),'None')
    poltable[cols[3]][badtask] = np.full(len(badtask),'None')
    for i in range(4,8):
        poltable[cols[i]][badtask] = np.full(len(badtask),1.0)
    for c in cols[8:]:
        poltable[c][badtask] = np.full(len(badtask),np.nan)

    print(poltable[badtask])
    
    ascii.write(poltable,
                os.path.join(filedir,'pol_allbeams.csv'),
                overwrite = True,format='csv')
    

#helper script to put updated continuum files in right place
def move_cont_valid():
    """
    Tom ran more continuum validation but provided extra files
    So I want to copy just the things I want to right location
    """
    valid_file_list = glob.glob(
       ("/Users/adams/Downloads/contQual/"
        "[1-2][0-9][0-1][0-9][0-3][0-9][0-9][0-9][0-9]/"
        "dynamicRange.dat"))
    contdir = os.path.join(filedir,"cont_valid")
    for task in valid_file_list:
        head, tail = os.path.split(task)
        tid = head[-9:]
        targetdir = os.path.join(contdir,tid)
        if not os.path.exists(targetdir):
            os.makedirs(targetdir)
        shutil.copy(task,targetdir)
