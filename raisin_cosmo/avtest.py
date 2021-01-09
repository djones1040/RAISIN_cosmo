#!/usr/bin/env python
# D. Jones - 12/2/20

import numpy as np
import pylab as plt
plt.ion()
from txtobj import txtobj
import cosmo

fr = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT000.FITRES',fitresheader=True)
frnb = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT000_nobiascor.FITRES',fitresheader=True)
frav1 = txtobj('output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',fitresheader=True)
frav2 = txtobj('output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT',fitresheader=True)
frav3 = txtobj('output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT',fitresheader=True)

def main():

    #plt.errorbar(fr.zCMB,fr.DLMAG-cosmo.mu(fr.zCMB),yerr=fr.DLMAGERR,fmt='o')
    #plt.errorbar(frnb.zCMB,frnb.DLMAG-cosmo.mu(frnb.zCMB),yerr=frnb.DLMAGERR,fmt='o')

    fr.optAV,fr.optAVERR = np.zeros(len(fr.CID)),np.zeros(len(fr.CID))
    for j,i in enumerate(fr.CID):
        if i in frav1.CID:
            fr.optAV[j] = frav1.AV[frav1.CID == i][0]
            fr.optAVERR[j] = frav1.AVERR[frav1.CID == i][0]
        elif i in frav2.CID:
            fr.optAV[j] = frav2.AV[frav2.CID == i][0]
            fr.optAVERR[j] = frav2.AVERR[frav2.CID == i][0]
        elif i in frav3.CID:
            fr.optAV[j] = frav3.AV[frav3.CID == i][0]
            fr.optAVERR[j] = frav3.AVERR[frav3.CID == i][0]
        else: raise RuntimeError

    frnb.optAV,frnb.optAVERR = np.zeros(len(frnb.CID)),np.zeros(len(frnb.CID))
    for j,i in enumerate(frnb.CID):
        if i in frav1.CID:
            frnb.optAV[j] = frav1.AV[frav1.CID == i][0]
            frnb.optAVERR[j] = frav1.AVERR[frav1.CID == i][0]
        elif i in frav2.CID:
            frnb.optAV[j] = frav2.AV[frav2.CID == i][0]
            frnb.optAVERR[j] = frav2.AVERR[frav2.CID == i][0]
        elif i in frav3.CID:
            frnb.optAV[j] = frav3.AV[frav3.CID == i][0]
            frnb.optAVERR[j] = frav3.AVERR[frav3.CID == i][0]
        else: raise RuntimeError

        
    plt.errorbar(fr.optAV,fr.DLMAG-cosmo.mu(fr.zCMB),xerr=fr.optAVERR,yerr=fr.DLMAGERR,fmt='o',label='with bias corr.')
    plt.errorbar(frnb.optAV,frnb.DLMAG-cosmo.mu(frnb.zCMB),xerr=frnb.optAVERR,yerr=frnb.DLMAGERR,fmt='o',label='without bias corr.')
    plt.xlabel('$A_V$',fontsize=15)
    plt.ylabel('Hubble Residual',fontsize=15)
    plt.legend()
    
    import pdb; pdb.set_trace()

if __name__ == "__main__":
    main()
