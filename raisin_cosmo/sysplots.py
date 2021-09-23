#!/usr/bin/env python
# D. Jones - 9/4/20

import pylab as plt
plt.ion()
import numpy as np
from txtobj import txtobj
from scipy.stats import binned_statistic
import cosmo

def main():
    ax1,ax2,ax3 = plt.subplot(221),plt.subplot(222),plt.subplot(223)
    ax4 = plt.subplot(224)

    lcbase = txtobj('output/cosmo_fitres/RAISIN_stat_lcparams.txt')
    lcpv = txtobj('output/cosmo_fitres/RAISIN_pecvel_lcparams.txt')
    lcpc = txtobj('output/cosmo_fitres/RAISIN_photcal_lcparams.txt')
    lcbc = txtobj('output/cosmo_fitres/RAISIN_biascor_lcparams.txt')
    lcms = txtobj('output/cosmo_fitres/RAISIN_massdivide_lcparams.txt')

    zbins = np.append(np.linspace(0.01,0.08,5),np.linspace(0.1,0.7,5))
    #pvbins = binned_statistic(lcbase.zcmb,lcpv.mb-lcbase.mb,bins=zbins,statistic='median').statistic
    #ax1.plot((zbins[1:]+zbins[:-1])/2.,pvbins,'o-',color='k')
    #ax1.set_title('Peculiar Velocities')

    #pcbins = binned_statistic(lcbase.zcmb,lcpc.mb-lcbase.mb,bins=zbins,statistic='median').statistic
    #ax2.plot((zbins[1:]+zbins[:-1])/2.,pcbins,'o-',color='k')
    #ax2.set_title('Phot. Cal.')

    #bcbins = binned_statistic(lcbase.zcmb,lcbc.mb-lcbase.mb,bins=zbins,statistic='median').statistic
    #ax3.plot((zbins[1:]+zbins[:-1])/2.,bcbins,'o-',color='k')
    #ax3.set_title('Bias Corr.')

    #msbins = binned_statistic(lcbase.zcmb,lcms.mb-lcbase.mb,bins=zbins,statistic='median').statistic
    #ax4.plot((zbins[1:]+zbins[:-1])/2.,msbins,'o-',color='k')
    #ax4.set_title('Mass Step')

    pvbins = binned_statistic(lcbase.zcmb,np.sqrt(lcpv.dmb**2.-lcbase.dmb**2.),bins=zbins,statistic='median').statistic
    ax1.plot((zbins[1:]+zbins[:-1])/2.,pvbins,'o-',color='k')
    ax1.set_title('Peculiar Velocities')

    pcbins = binned_statistic(lcbase.zcmb,np.sqrt(lcpc.dmb**2.-lcbase.dmb**2.),bins=zbins,statistic='median').statistic
    ax2.plot((zbins[1:]+zbins[:-1])/2.,pcbins,'o-',color='k')
    ax2.set_title('Phot. Cal.')

    bcbins = binned_statistic(lcbase.zcmb,np.sqrt(lcbc.dmb**2.-lcbase.dmb**2.),bins=zbins,statistic='median').statistic
    ax3.plot((zbins[1:]+zbins[:-1])/2.,bcbins,'o-',color='k')
    ax3.set_title('Bias Corr.')

    msbins = binned_statistic(lcbase.zcmb,np.sqrt(lcms.dmb**2.-lcbase.dmb**2.),bins=zbins,statistic='median').statistic
    ax4.plot((zbins[1:]+zbins[:-1])/2.,msbins,'o-',color='k')
    ax4.set_title('Mass Step')
    
    
    import pdb; pdb.set_trace()

def main_paper():

    plt.subplots_adjust(hspace=0,wspace=0,left=0.2,right=0.97)
    
    ax1,ax2,ax3 = plt.subplot(331),plt.subplot(332),plt.subplot(333)
    ax4,ax5,ax6 = plt.subplot(334),plt.subplot(335),plt.subplot(336)
    ax7,ax8,ax9 = plt.subplot(337),plt.subplot(338),plt.subplot(339)
    for ax in [ax7,ax8,ax9]:
        ax.set_xlabel('$z_{CMB}$',fontsize=15)
    for ax in [ax1,ax4,ax7]:
        ax.set_ylabel('$\delta\mu$ (mag)',fontsize=15)
    for ax in [ax2,ax3,ax5,ax6]:
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
    ax1.xaxis.set_ticklabels([])
    ax4.xaxis.set_ticklabels([])
    ax8.yaxis.set_ticklabels([])
    ax9.yaxis.set_ticklabels([])
    
    for ax in [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]:
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        ax.set_xlim([0.0,0.6])
        ax.set_ylim([-0.027,0.027])
        ax.axhline(0.0,lw=2,color='k')
    for ax in [ax7,ax8,ax9]:
        ax.set_xlim([0.0,0.6])
        #ax.set_ylim([-0.1,0.1])
        ax.set_ylim([-0.06,0.06])
        ax.yaxis.set_ticks([-0.04,-0.02,0.0,0.02,0.04])
    # plots
    # cal:     HST_CAL CSP_Y CSP_J
    # biascor: biascor_av_lowz biascor_av_highz biascor_shape_highz
    # other:   mass_step vpec kcor

    lcbase = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT000_new.FITRES',fitresheader=True)
    zbins = np.linspace(0.01,0.8,15)
    zplot = (zbins[1:]+zbins[:-1])/2.
    histnum = np.histogram(lcbase.zHD,bins=zbins)[0]
    zbins = np.append(zbins[:-1][histnum > 2],zbins[-1])
    zplot = zplot[histnum > 2]
    
    # cal: FITOPT002 FITOPT006 FITOPT007
    for ax,name in zip([ax1,ax2,ax3],['$HST$','CSP $Y$','CSP $J$']):
        ax.text(0.5,0.9,name,transform=ax.transAxes,ha='center',va='center')
    plt.gcf().text(0.03,0.77,'Cal.',ha='center',va='center',rotation=90,fontsize=20)
    plt.gcf().text(0.03,0.5,'Bias Corr.',ha='center',va='center',rotation=90,fontsize=20)
    plt.gcf().text(0.03,0.23,'Misc.',ha='center',va='center',rotation=90,fontsize=20)
    
    frhst = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT002_new.FITRES',fitresheader=True)
    frY = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT007_new.FITRES',fitresheader=True)
    frJ = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT008_new.FITRES',fitresheader=True)

    def sys_average(x,fr=None,frbase=None):
        global_off = np.average(fr.DLMAG-frbase.DLMAG,weights=1/(fr.DLMAGERR**2.))
        dm2=(fr.DLMAG[x]-cosmo.mu(fr.zHD[x]))-global_off-(frbase.DLMAG[x]-cosmo.mu(frbase.zHD[x]))
        return np.average(dm2,weights=1/fr.DLMAGERR[x]**2.)

    hstbins = binned_statistic(lcbase.zHD,range(len(lcbase.zHD)),bins=zbins,statistic=lambda values: sys_average(values,frhst,lcbase)).statistic
    ax1.plot(zplot[hstbins == hstbins],hstbins[hstbins == hstbins],'o-',color='0.6')
    Ybins = binned_statistic(lcbase.zHD,range(len(lcbase.zHD)),bins=zbins,statistic=lambda values: sys_average(values,frY,lcbase)).statistic
    ax2.plot(zplot[Ybins == Ybins],Ybins[Ybins == Ybins],'o-',color='0.6')
    Jbins = binned_statistic(lcbase.zHD,range(len(lcbase.zHD)),bins=zbins,statistic=lambda values: sys_average(values,frJ,lcbase)).statistic
    ax3.plot(zplot[Jbins == Jbins],Jbins[Jbins == Jbins],'o-',color='0.6')

    # biascor
    #FITOPT009 FITOPT011 FITOPT012
    frb1 = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT011_new.FITRES',fitresheader=True)
    frb2 = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT013_new.FITRES',fitresheader=True)
    frb3 = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT012_new.FITRES',fitresheader=True)
    b1bins = binned_statistic(lcbase.zHD,range(len(lcbase.zHD)),bins=zbins,statistic=lambda values: sys_average(values,frb1,lcbase)).statistic
    ax4.plot(zplot[b1bins == b1bins],b1bins[b1bins == b1bins],'o-',color='0.6')
    b2bins = binned_statistic(lcbase.zHD,range(len(lcbase.zHD)),bins=zbins,statistic=lambda values: sys_average(values,frb2,lcbase)).statistic
    ax5.plot(zplot[b2bins == b2bins],b2bins[b2bins == b2bins],'o-',color='0.6')
    b3bins = binned_statistic(lcbase.zHD,range(len(lcbase.zHD)),bins=zbins,statistic=lambda values: sys_average(values,frb3,lcbase)).statistic
    ax6.plot(zplot[b3bins == b3bins],b3bins[b3bins == b3bins],'o-',color='0.6')
    for ax,name in zip([ax4,ax5,ax6],['CSP, 1$\sigma$ $A_V$','high-$z$, 1$\sigma$ $A_V$','high-$z$, 1$\sigma$ $s_{BV}$']):
        ax.text(0.5,0.9,name,transform=ax.transAxes,ha='center',va='center')

    # other
    # FITOPT005 FITOPT003 FITOPTxxx?
    frb1 = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT005_new.FITRES',fitresheader=True)
    frb2 = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT006_new.FITRES',fitresheader=True)
    frb3 = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT014_new.FITRES',fitresheader=True)
    b1bins = binned_statistic(lcbase.zHD,range(len(lcbase.zHD)),bins=zbins,statistic=lambda values: sys_average(values,frb1,lcbase)).statistic
    ax7.plot(zplot[b1bins == b1bins],b1bins[b1bins == b1bins],'o-',color='0.6')
    b2bins = binned_statistic(lcbase.zHD,range(len(lcbase.zHD)),bins=zbins,statistic=lambda values: sys_average(values,frb2,lcbase)).statistic
    ax8.plot(zplot[b2bins == b2bins],b2bins[b2bins == b2bins],'o-',color='0.6')
    b3bins = binned_statistic(lcbase.zHD,range(len(lcbase.zHD)),bins=zbins,statistic=lambda values: sys_average(values,frb3,lcbase)).statistic
    ax9.plot(zplot[b3bins == b3bins],b3bins[b3bins == b3bins],'o-',color='0.6')
    for ax,name in zip([ax7,ax8,ax9],['1$\sigma$ mass step','NIR SN model','$K$-corr.']):
        ax.text(0.5,0.9,name,transform=ax.transAxes,ha='center',va='center',bbox={'facecolor':'1.0','edgecolor':'1.0','pad':0.0,'alpha':0.7})
    
    import pdb; pdb.set_trace()

    
if __name__ == "__main__":
    #main()
    main_paper()
