#!/usr/bin/env python

import pylab as plt
plt.ion()
import numpy as np

def main():
    plt.rcParams['figure.figsize'] = (8,8)
    plt.subplots_adjust(hspace=1.0,top=0.93,bottom=0.07)#,wspace=0)
    
    specfiles = ['ps1-450339.dat','ps1-480464.dat',
                 'ps1-490037.dat','ps1-490521.dat',
                 'ps1-520188.dat',
                 'des15x2nkz.dat','des16c2cva.dat',
                 'des16s2afz.dat']
    snidfiles = ['ps1-450339_snidflux.dat','ps1-480464_snidflux.dat',
                 'ps1-490037_snidflux.dat','ps1-490521_snidflux.dat',
                 'ps1-520188_snidflux.dat',
                 'des15x2nkz_snidflux.dat','des16c2cva_snidflux.dat',
                 'des16s2afz_snidflux.dat']
    snidlabels = ['SN 2003du (+0 days)','SN 2007F (+3 days)','SN 2005na (-2 days)',
              'SN 2006ot (+3 days)','SN 1999dq (-9 days)','SN 2005cf (-1 day)',
              'SN 2004L (+3 days)','SN 2006S (-6 days)']
    zlist = [0.4096,0.2212,0.4231,0.3412,0.2804,0.468800,0.402900,0.483000]
    for sf,snf,sl,i,z in zip(specfiles,snidfiles,snidlabels,range(len(specfiles)),zlist):

        ax = plt.subplot(f"42{i+1:.0f}")
        ax.tick_params(top="off",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        
        wave,flux = np.loadtxt(sf,unpack=True)
        swave,sflux = np.loadtxt(snf,unpack=True)

        ax.plot(wave,flux,color='k',
                label=sf.split('_')[0].replace('ps1','PS1').replace('des','DES').replace('.dat',''))
        ax.plot(swave,sflux,color='r',
                label=sl)
        ax.legend(loc='upper right',prop={'size':8})
        #if i >= 6:
        ax.set_xlabel(r'Observed wavelength (${\rm \AA}$)',labelpad=0)
        #if i in [0,2,4,6]:
        ax.set_ylabel('Normalized Flux')
        #ax.set_xlim([4000,9000])
        ax.set_ylim([0.3,2.3])

        axtop = ax.twiny()
        axtop.xaxis.tick_top()
        axtop.xaxis.set_ticks([4000*(1+z),5000*(1+z),6000*(1+z),7000*(1+z)])
        axtop.xaxis.set_ticklabels([4000,5000,6000,7000])
        axtop.set_xlim(ax.get_xlim())
        axtop.tick_params(direction="inout",length=8, width=1.5)
        axtop.set_xlabel(r'Rest Wavelength (${\rm \AA}$)')
        #ax.set_xticks([5000,6000,7000,8000])
        
    import pdb; pdb.set_trace()
            
if __name__ == "__main__":
    main()
