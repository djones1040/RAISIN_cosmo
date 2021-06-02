#!/usr/bin/env python

import pylab as plt
plt.ion()
import numpy as np

def main():
    plt.rcParams['figure.figsize'] = (8,8)
    plt.subplots_adjust(hspace=0.5)#,wspace=0)
    
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
    for sf,snf,sl,i in zip(specfiles,snidfiles,snidlabels,range(len(specfiles))):

        ax = plt.subplot(f"42{i+1:.0f}")
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        
        wave,flux = np.loadtxt(sf,unpack=True)
        swave,sflux = np.loadtxt(snf,unpack=True)

        ax.plot(wave,flux,color='k',
                label=sf.split('_')[0].replace('ps1','PS1').replace('des','DES'))
        ax.plot(swave,sflux,color='r',
                label=sl)
        ax.legend(loc='upper right')
        #if i >= 6:
        ax.set_xlabel(r'Observed wavelength (${\rm \AA}$)')
        #if i in [0,2,4,6]:
        ax.set_ylabel('Normalized Flux')
        #ax.set_xlim([4000,9000])
        ax.set_ylim([0.3,2.3])
        #ax.set_xticks([5000,6000,7000,8000])
        
    import pdb; pdb.set_trace()
            
if __name__ == "__main__":
    main()
