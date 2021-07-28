#!/usr/bin/env python
# D. Jones - 6/29/21

import numpy as np
import pylab as plt
plt.ion()
from txtobj import txtobj
import astropy.table as at

import palettable

from palettable.colorbrewer.qualitative import Dark2_8 as palettable_color

def lowz():
    # host masses from lowz
    fr = txtobj('output/fitres_cosmo/CSP.FITRES',fitresheader=True)
    fr = txtobj('hosts/lowz_optical_test.txt')
    
    # host masses from Roman
    data = at.Table.read('hosts/snprop_roman.dat',format='ascii',data_start=2)

    #hostmass_dj,hostmass_r,snid = [],[],[]
    #for j,i in enumerate(fr.CID):
    #    if 'sn'+i in data['Name']:
    #        hostmass_dj += [fr.HOST_LOGMASS[j]]
    #        hostmass_r += [data['logMs'][data['Name'] == 'sn'+i][0]]
    #        snid += [i]

    hostmass_dj,hostmass_r,snid = [],[],[]
    for j,i in enumerate(fr.snid):
        if 'sn'+i in data['Name']:
            hostmass_dj += [fr.logmass[j]]
            hostmass_r += [data['logMs'][data['Name'] == 'sn'+i][0]]
            snid += [i]

            
    hostmass_dj,hostmass_r = np.array(hostmass_dj),np.array(hostmass_r)
    snid = np.array(snid)
    
    plt.plot(hostmass_dj,hostmass_r,'o')
    plt.plot([9,12],[9,12],color='k')
    
    import pdb; pdb.set_trace()

def des():
    # host masses from lowz
    fr = txtobj('output/fitres_cosmo/DES.FITRES',fitresheader=True)

    # host masses from Roman
    data = at.Table.read('hosts/RAISINS_Masses_des.csv',format='ascii')

    hostmass_dj,hostmass_des,snid = [],[],[]
    for j,i in enumerate(fr.CID):
        if i in data['TRANSIENT_NAME']:
            hostmass_dj += [fr.HOST_LOGMASS[j]]
            hostmass_des += [data['MASS'][data['TRANSIENT_NAME'] == i][0]]
            snid += [i]
            
    hostmass_dj,hostmass_des = np.array(hostmass_dj),np.array(hostmass_des)
    snid = np.array(snid)
    
    plt.plot(hostmass_dj,hostmass_des,'o')
    plt.plot([8,12],[8,12],color='k')
    
    import pdb; pdb.set_trace()

def histplot():
    fr = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT000.FITRES',fitresheader=True)

    ax = plt.axes()
    ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    ax.set_ylabel(r'N$_{\rm SNe}$',fontsize=15)
    ax.set_xlabel(r'log($M_{\ast}/M_{\odot}$)',fontsize=15,labelpad=0)

    massbins = np.linspace(7,13,12)
    ax.hist(fr.HOST_LOGMASS[fr.IDSURVEY == 5],bins=massbins,ec='k',label='CSP',alpha=0.5,hatch='\\')
    ax.hist(fr.HOST_LOGMASS[fr.IDSURVEY == 15],bins=massbins,ec='k',label='PS1 (RAISIN1)',alpha=0.5)
    ax.hist(fr.HOST_LOGMASS[fr.IDSURVEY == 10],bins=massbins,ec='k',label='DES (RAISIN2)',alpha=0.5,hatch='//')
    ax.set_xlim([7.5,12])
    ax.yaxis.set_ticks([0,5,10,15])
    
    ax.legend(prop={'size':12})
    
    import pdb; pdb.set_trace()

def cdf():
    fr = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT000.FITRES',fitresheader=True)

    ax = plt.axes()
    ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    ax.set_ylabel(r'Cumulative Fraction',fontsize=15)
    ax.set_xlabel(r'log($M_{\ast}/M_{\odot}$)',fontsize=15,labelpad=0)

    massbins = np.linspace(7,13,50)
    lowz_hist = np.histogram(fr.HOST_LOGMASS[fr.IDSURVEY == 5],bins=massbins) #,ec='k',label='CSP',alpha=0.5,hatch='\\')
    ps1_hist = np.histogram(fr.HOST_LOGMASS[fr.IDSURVEY == 15],bins=massbins) #,ec='k',label='PS1 (RAISIN1)',alpha=0.5)
    des_hist = np.histogram(fr.HOST_LOGMASS[fr.IDSURVEY == 10],bins=massbins) #,ec='k',label='DES (RAISIN2)',alpha=0.5,hatch='//')
    ax.plot(lowz_hist[1][:-1],
			 np.cumsum(lowz_hist[0])/float(np.sum(lowz_hist[0])),drawstyle='steps-mid',
			 label='CSP',color=palettable_color.mpl_colors[0],lw=2)
    ax.plot(lowz_hist[1][:-1],
			 np.cumsum(ps1_hist[0])/float(np.sum(ps1_hist[0])),drawstyle='steps-mid',
			 label='PS1 (RAISIN1)',color=palettable_color.mpl_colors[1],lw=2)
    ax.plot(lowz_hist[1][:-1],
			 np.cumsum(des_hist[0])/float(np.sum(des_hist[0])),drawstyle='steps-mid',
			 label='DES (RAISIN2)',color=palettable_color.mpl_colors[2],lw=2)
    
    ax.set_xlim([8,11.5])
    #ax.yaxis.set_ticks([0,5,10,15])
    
    ax.legend(prop={'size':12})
    ax.set_ylim([0,1.05])
    import pdb; pdb.set_trace()

    
if __name__ == "__main__":
    #lowz()
    #des()
    #histplot()
    cdf()
