#!/usr/bin/env python
# D. Jones - 12/17/20

import numpy as np
import pylab as plt
from matplotlib import colors
from astropy.io import fits
plt.ion()
import snana
import glob

def main():

    plt.subplots_adjust(hspace=0,wspace=0,right=0.8)
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)
    for ax in [ax1,ax2,ax3,ax4]:
        ax.tick_params(top="on",bottom="on",left="on",right="on",
                       direction="inout",length=8, width=1.5)
    cax = plt.axes([0.85,0.11,0.05,0.77])
    
    snpyfits = fits.open('snoopy.B18/SNooPy_test.fits')
    filters = snpyfits[0].header['FILTERS']
    
    NGRID_FILT = snpyfits[1].data['NBIN'][4]
    NGRID_TREST = snpyfits[1].data['NBIN'][5]
    NWDPAD = 4
    ST_GRID = np.linspace(0.7,1.3,snpyfits[1].data['NBIN'][snpyfits[1].data['PAR_NAME'] == 'STRETCH'][0])

    I2MAG = snpyfits[9].data['I2MAG']
    I2MAGERR = snpyfits[9].data['I2MAGERR']
    TREST = snpyfits[7].data['TREST']
    
    ILCOFF = [1, 1, 1, 1]
    norm=colors.Normalize(vmin=0.7,vmax=1.3)
    cmap=plt.get_cmap('RdBu')
    for st in ST_GRID:
        
        # redshift, AV, RV, luminosity
        INDX = [1,1,1,np.where(ST_GRID == st)[0]+1]
        ILC = 1
        for i in range(4):
            ILC += ILCOFF[i] * (INDX[i] - 1)

        PTR = 1 + ((NGRID_FILT * NGRID_TREST) + NWDPAD) * (ILC-1)
        
        for iflt,flt in enumerate(filters):
            
            lctmpl = snpyfits[9].data['I2MAG'][int(PTR)+1+NGRID_TREST*iflt:int(PTR)+1+NGRID_TREST*(iflt+1)]/1000.
            color=cmap(norm(st))
            if flt == 'i': ax1.plot(TREST,lctmpl,color=color)
            if flt == 'Y': ax2.plot(TREST,lctmpl,color=color)
            if flt == 'J': ax3.plot(TREST,lctmpl,color=color)
            if flt == 'H': ax4.plot(TREST,lctmpl,color=color)
            
            #snpyfits[9].data['I2MAGERR'][int(PTR)+1+NGRID_TREST*iflt:int(PTR)+1+NGRID_TREST*(iflt+1)] = \
            #   interptmplerr*1000.
            #if flt in 'iYJH' and st > 1.005 and st < 1.01:
            #    print(st,flt,np.min(lctmpl))
            #    import pdb; pdb.set_trace()
            #if flt in 'iYJH' and st > 1.06 and st < 1.07:
            #    print(st,flt,np.min(lctmpl))
            #    import pdb; pdb.set_trace()
            #if flt in 'iYJH' and st > 1.075 and st < 1.08:

                #import pdb; pdb.set_trace()
            #if flt == 'iYJH' and st > 1.195 and st < 1.2:
            #    print(st,flt,np.min(lctmpl))
                #import pdb; pdb.set_trace()
            #if flt == 'i' and st > 1.175 and st < 1.18: import pdb; pdb.set_trace()
            #if flt == 'i' and st > 1.295 and st < 1.3: import pdb; pdb.set_trace()
            
    for ax,flt in zip([ax1,ax2,ax3,ax4],'iYJH'):
        #ax.invert_yaxis()
        ax.text(0.5,0.9,f"${flt}$",va='center',ha='center',transform=ax.transAxes,fontsize=15)
        ax.set_ylim([-16.5,-19.2])
        ax.set_xlim([-16,50])
        ax.yaxis.set_ticks([-17,-18,-19])
        
    #ax3.set_xlabel('Phase (days)')
    xxl = ax4.set_xlabel('Phase (days)',fontsize=15)
    xxl.set_position((xxl.get_position()[1],1))
    xxl.set_horizontalalignment('center')
    
    #ax1.set_ylabel('Absolute Mag')
    #ax3.set_ylabel('Absolute Mag')
    yyl = ax3.set_ylabel('Absolute Mag',fontsize=15,labelpad=10)
    yyl.set_position((yyl.get_position()[0],1))
    yyl.set_verticalalignment('center')
    ax2.yaxis.set_ticklabels([])
    ax4.yaxis.set_ticklabels([])
    
    sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap)
    plt.gcf().colorbar(sm,cax=cax)
    cax.set_ylabel('$s_{BV}$',fontsize=15)


    goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                               np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                               np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)))
    raisin_lcs = np.append(glob.glob('data/Photometry/PS1_RAISIN/*.snana.dat'),
                           glob.glob('data/Photometry/DES_RAISIN/*.snana.dat'))
    lowz_lcs = glob.glob('data/Photometry/CSPDR3_RAISIN/*.PKMJD.DAT')

    raisin_phases = []
    for r in raisin_lcs:
        sn = snana.SuperNova(r)
        if sn.SNID in goodcids:
            phases = (sn.MJD[(sn.FLT == 'H')]-sn.PEAKMJD)/(1+float(sn.REDSHIFT_FINAL.split()[0]))
            raisin_phases = np.append(raisin_phases,phases)
    raisin_min_phase = np.sort(raisin_phases)[int(len(raisin_phases)*0.16)]
    raisin_max_phase = np.sort(raisin_phases)[int(len(raisin_phases)*0.84)]

    lowz_phases = []
    for r in lowz_lcs:
        sn = snana.SuperNova(r)
        if sn.SNID in goodcids:
            phases = (sn.MJD[(sn.FLT == 'H') | (sn.FLT == 'h') | (sn.FLT == 'J') | (sn.FLT == 'j') | (sn.FLT == 'Y') | (sn.FLT == 'y')]-sn.PEAKMJD)/(1+float(sn.REDSHIFT_FINAL.split()[0]))
            lowz_phases = np.append(lowz_phases,phases[phases < 45])
    lowz_min_phase = np.sort(lowz_phases)[int(len(lowz_phases)*0.16)]
    lowz_max_phase = np.sort(lowz_phases)[int(len(lowz_phases)*0.84)]

    
    # RAISIN phases, CSP phases
    for ax in [ax1,ax2,ax3,ax4]:
        ax.axvspan(
            raisin_min_phase,raisin_max_phase,
            color='0.5',alpha=0.5)
        #ax.axvspan(
        #    lowz_min_phase,lowz_max_phase,
        #    color='green',alpha=0.5)
    
    import pdb; pdb.set_trace()

def main_flux():

    plt.subplots_adjust(hspace=0,wspace=0,right=0.8)
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)
    for ax in [ax1,ax2,ax3,ax4]:
        ax.tick_params(top="on",bottom="on",left="on",right="on",
                       direction="inout",length=8, width=1.5)
    cax = plt.axes([0.85,0.11,0.05,0.77])
    
    snpyfits = fits.open('snoopy.B18/SNooPy_test.fits')
    filters = snpyfits[0].header['FILTERS']
    
    NGRID_FILT = snpyfits[1].data['NBIN'][4]
    NGRID_TREST = snpyfits[1].data['NBIN'][5]
    NWDPAD = 4
    ST_GRID = np.linspace(0.7,1.3,snpyfits[1].data['NBIN'][snpyfits[1].data['PAR_NAME'] == 'STRETCH'][0])

    I2MAG = snpyfits[9].data['I2MAG']
    I2MAGERR = snpyfits[9].data['I2MAGERR']
    TREST = snpyfits[7].data['TREST']
    
    ILCOFF = [1, 1, 1, 1]
    norm=colors.Normalize(vmin=0.7,vmax=1.3)
    cmap=plt.get_cmap('RdBu')
    for st in ST_GRID:
        
        # redshift, AV, RV, luminosity
        INDX = [1,1,1,np.where(ST_GRID == st)[0]+1]
        ILC = 1
        for i in range(4):
            ILC += ILCOFF[i] * (INDX[i] - 1)

        PTR = 1 + ((NGRID_FILT * NGRID_TREST) + NWDPAD) * (ILC-1)
        
        for iflt,flt in enumerate(filters):
            
            lctmpl = snpyfits[9].data['I2MAG'][int(PTR)+1+NGRID_TREST*iflt:int(PTR)+1+NGRID_TREST*(iflt+1)]/1000.
            color=cmap(norm(st))
            if flt == 'i': ax1.plot(TREST,10**(-0.4*(lctmpl-27.5)),color=color)
            if flt == 'Y': ax2.plot(TREST,10**(-0.4*(lctmpl-27.5)),color=color)
            if flt == 'J': ax3.plot(TREST,10**(-0.4*(lctmpl-27.5)),color=color)
            if flt == 'H': ax4.plot(TREST,10**(-0.4*(lctmpl-27.5)),color=color)
            
            
    for ax,flt in zip([ax1,ax2,ax3,ax4],'iYJH'):
        #ax.invert_yaxis()
        ax.text(0.5,0.9,f"${flt}$",va='center',ha='center',transform=ax.transAxes,fontsize=15)
        #ax.set_ylim([-16.5,-19.2])
        #ax.set_xlim([-16,50])
        #ax.yaxis.set_ticks([-17,-18,-19])
        
    ax3.set_xlabel('Phase')
    ax4.set_xlabel('Phase')
    ax1.set_ylabel('Mag')
    ax3.set_ylabel('Mag')
    ax2.yaxis.set_ticklabels([])
    ax4.yaxis.set_ticklabels([])
    
    sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap)
    plt.gcf().colorbar(sm,cax=cax)
    cax.set_ylabel('$s_{BV}$',fontsize=15)
    
    import pdb; pdb.set_trace()

    
if __name__ == "__main__":
    main()
    #main_flux()
