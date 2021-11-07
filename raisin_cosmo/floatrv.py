#!/usr/bin/env python
# D. Jones - 10/21/21

import numpy as np
from txtobj import txtobj
import pylab as plt
plt.ion()
import cosmo

_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',dtype=str),
                            np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',dtype=str),
                            np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',dtype=str)))

def apply_all_cuts(fr,fropt,restrict_to_good_list=True):

    # AV
    iGoodAV = np.zeros(len(fr.CID),dtype=bool)
    for j,i in enumerate(fr.CID):
        if i in fropt.CID and fropt.AV[fropt.CID == i] < 0.3*fropt.RV[fropt.CID == i]:
            iGoodAV[j] = True

    # reasonable stretch
    iGoodSt = np.zeros(len(fr.CID),dtype=bool)
    for j,i in enumerate(fr.CID):
        if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > 0.75 and fropt.STRETCH[fropt.CID == i] < 1.185 and \
           fropt.STRETCHERR[fropt.CID == i] < 0.3: #175:
            iGoodSt[j] = True

    #iGoodSt = (fropt.STRETCH > 0.8) & (fropt.STRETCH < 1.3)

    #for k in fr.__dict__.keys():
    #    fr.__dict__[k] = fr.__dict__[k][iGoodAV & iGoodSt]

    #iGoodRV = (fropt.RVERR < 0.5) & (fropt.RV > 0)

    #for k in fr.__dict__.keys():
    #    fr.__dict__[k] = fr.__dict__[k][iGoodRV]

    
    #restrict_to_good_list = False
    if restrict_to_good_list:
        # then case-by-case removal of bad stuff
        # because some RAISIN things have bad photometry
        iGood = np.array([],dtype=int)
        for j,i in enumerate(fr.CID):
            if i in _goodcids: iGood = np.append(iGood,j)
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]

    return fr

def main():
    fr1 = txtobj('output/fit_optical/CSP_RAISIN_optnir_FLOATRV.FITRES.TEXT',fitresheader=True)
    fr2 = txtobj('output/fit_optical/PS1_RAISIN_optnir_FLOATRV.FITRES.TEXT',fitresheader=True)
    fr3 = txtobj('output/fit_optical/DES_RAISIN_optnir_FLOATRV.FITRES.TEXT',fitresheader=True)

    fr1f = txtobj('output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',fitresheader=True)
    fr2f = txtobj('output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT',fitresheader=True)
    fr3f = txtobj('output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT',fitresheader=True)

    fr1.RV = fr1.AV/fr1.AVRV
    fr2.RV = fr2.AV/fr2.AVRV
    fr3.RV = fr3.AV/fr3.AVRV
    
    fr1 = apply_all_cuts(fr1,fr1)
    fr2 = apply_all_cuts(fr2,fr2)
    fr3 = apply_all_cuts(fr3,fr3)
    fr1f = apply_all_cuts(fr1f,fr1f)
    fr2f = apply_all_cuts(fr2f,fr2f)
    fr3f = apply_all_cuts(fr3f,fr3f)

    fr1.resid = fr1.DLMAG - cosmo.mu(fr1.zHD)
    fr2.resid = fr2.DLMAG - cosmo.mu(fr2.zHD)
    fr3.resid = fr3.DLMAG - cosmo.mu(fr3.zHD)
    fr1f.resid = fr1f.DLMAG - cosmo.mu(fr1f.zHD)
    fr2f.resid = fr2f.DLMAG - cosmo.mu(fr2f.zHD)
    fr3f.resid = fr3f.DLMAG - cosmo.mu(fr3f.zHD)
    
    #plt.hist(fr1.RV,histtype='step')
    #plt.hist(fr2.RV,histtype='step')
    #plt.hist(fr3.RV,histtype='step')
    
    iGood = fr1.AVERR/fr1.RV < 0.5
    plt.errorbar(fr1.AV[iGood]/fr1.RV[iGood],fr1.RV[iGood],
#                 xerr=fr1.RVERR[iGood],
                 yerr=fr1.AVERR[iGood]/fr1.RV[iGood],fmt='o',label='CSP')
    iGood = fr2.AVERR/fr2.RV < 0.5
    plt.errorbar(fr2.AV[iGood]/fr2.RV[iGood],fr2.RV[iGood],
#                 xerr=fr2.RVERR[iGood],
                 yerr=fr2.AVERR[iGood]/fr2.RV[iGood],fmt='o',label='RAISIN1')
    iGood = fr3.AVERR/fr3.RV < 0.5
    plt.errorbar(fr3.AV[iGood]/fr3.RV[iGood],fr3.RV[iGood],
#                 xerr=fr3.RVERR[iGood],
                 yerr=fr3.AVERR[iGood]/fr3.RV[iGood],fmt='o',label='RAISIN2')

    plt.xlabel('$E(B-V)$')
    plt.ylabel('$R_V$')
    plt.legend()
    import pdb; pdb.set_trace()
    
if __name__ == "__main__":
    main()
