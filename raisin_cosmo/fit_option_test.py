#!/usr/bin/env python
# D. Jones - 2/17/21

import pylab as plt
plt.ion()
import numpy as np
import cosmo
from txtobj import txtobj
from astropy.stats import sigma_clipped_stats
import copy

goodcids_raisin1 = np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',dtype=str,unpack=True)
goodcids_raisin2 = np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',dtype=str,unpack=True)
goodcids_lowz_csp = np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',dtype=str,unpack=True)

def mkcuts(fr,fr2,frlowz):
    # get the distance moduli
    iGood = np.array([],dtype=int)
    for j,i in enumerate(fr.CID):
        if i in goodcids_raisin1: iGood = np.append(iGood,j)
    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iGood]

    iGood = np.array([],dtype=int)
    for j,i in enumerate(fr2.CID):
        if i in goodcids_raisin2: iGood = np.append(iGood,j)
    for k in fr2.__dict__.keys():
        fr2.__dict__[k] = fr2.__dict__[k][iGood]

    iGood = np.array([],dtype=int)
    for j,i in enumerate(frlowz.CID):
        if i in goodcids_lowz_csp: iGood = np.append(iGood,j)
    for k in frlowz.__dict__.keys():
        frlowz.__dict__[k] = frlowz.__dict__[k][iGood]

    return fr,fr2,frlowz

def main():
    plt.rcParams['figure.figsize'] = (5,5)
    residbins = np.linspace(-1,1,15)
    plt.subplots_adjust(hspace=0,wspace=0,right=0.97,top=0.97)

    
    for snanafile_raisin1,snanafile_raisin2,snanafile_lowz_csp,variant,ax,i in \
        zip(['output/fit_nir/PS1_RAISIN.FITRES.TEXT',
             'output/fit_nir/PS1_RAISIN_NIR_COLOR.FITRES.TEXT',
             'output/fit_nir/PS1_RAISIN_NIR_SHAPE.FITRES.TEXT',
             'output/fit_nir/PS1_RAISIN_NIR_SHAPECOLOR.FITRES.TEXT'][::-1],
            ['output/fit_nir/DES_RAISIN.FITRES.TEXT',
             'output/fit_nir/DES_RAISIN_NIR_SHAPE.FITRES.TEXT',
             'output/fit_nir/DES_RAISIN_NIR_COLOR.FITRES.TEXT',
             'output/fit_nir/DES_RAISIN_NIR_SHAPECOLOR.FITRES.TEXT'][::-1],
            ['output/fit_nir/DES_RAISIN.FITRES.TEXT',
             'output/fit_nir/CSP_RAISIN_NIR_SHAPE.FITRES.TEXT',
             'output/fit_nir/CSP_RAISIN_NIR_COLOR.FITRES.TEXT',
             'output/fit_nir/CSP_RAISIN_NIR_SHAPECOLOR.FITRES.TEXT'][::-1],
            ['all fixed','$s_{{BV}}$','$A_V$','$s_{B{V}}$ and $A_V$'][::-1],
            [plt.subplot(411),plt.subplot(412),plt.subplot(413),plt.subplot(414)],
            [0,1,2,3]):

        #if i < 2: continue
        
        fr = txtobj(snanafile_raisin1,fitresheader=True)
        fr2 = txtobj(snanafile_raisin2,fitresheader=True)
        frlowz = txtobj(snanafile_lowz_csp,fitresheader=True)
        fr,fr2,frlowz = mkcuts(fr,fr2,frlowz)

        fr.resid = fr.DLMAG - cosmo.mu(fr.zHD)
        fr2.resid = fr2.DLMAG - cosmo.mu(fr2.zHD)
        frlowz.resid = frlowz.DLMAG - cosmo.mu(frlowz.zHD)
        #median_resid = np.median(np.concatenate((fr.resid,fr2.resid,frlowz.resid)))
        fr.resid -= np.median(fr.resid); fr2.resid -= np.median(fr2.resid); frlowz.resid -= np.median(frlowz.resid)

        resid_highz = np.append(fr.resid,fr2.resid)
        stretch_highz = np.append(fr.STRETCH,fr2.STRETCH)
        residerr_highz = np.append(fr.DLMAGERR,fr2.DLMAGERR)
        stretcherr_highz = np.append(fr.STRETCHERR,fr2.STRETCHERR)


        fr3 = txtobj('output/fit_nir/PS1_RAISIN.FITRES.TEXT',fitresheader=True)
        fr4 = txtobj('output/fit_nir/DES_RAISIN.FITRES.TEXT',fitresheader=True)
        idx3 = np.array([],dtype=int)
        for j in fr.CID:
            idx3 = np.append(idx3,np.where(fr3.CID == j)[0])
        for k in fr3.__dict__.keys():
            fr3.__dict__[k] = fr3.__dict__[k][idx3]
        idx4 = np.array([],dtype=int)
        for j in fr2.CID:
            idx4 = np.append(idx4,np.where(fr4.CID == j)[0])
        for k in fr4.__dict__.keys():
            fr4.__dict__[k] = fr4.__dict__[k][idx4]

        #resid_highz = np.append(fr3.DLMAG - cosmo.mu(fr3.zHD),fr4.DLMAG - cosmo.mu(fr4.zHD))
        #residerr_highz = np.append(fr3.DLMAGERR,fr4.DLMAGERR)
        #plt.errorbar(stretch_highz[stretcherr_highz < 0.2],resid_highz[stretcherr_highz < 0.2],xerr=stretcherr_highz[stretcherr_highz < 0.2],yerr=residerr_highz[stretcherr_highz < 0.2],fmt='.')

        
        ax.hist(resid_highz,histtype='stepfilled',
                label=f'{variant}, {len(resid_highz):.0f} SNe\nRMS = {np.std(resid_highz):.3f} mag',
                alpha=1.0,bins=residbins,lw=2,color=f'C{i}')
        #ax = plt.axes()
        ax.set_ylabel(r'N$_{\rm SNe}$',fontsize=15)
        ax.set_ylim([0,15])
        ax.yaxis.set_ticks([5,10])
        if i == 3: ax.set_xlabel('Hubble Residual')
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        
        ax.legend(loc='upper left')
        #import pdb; pdb.set_trace()
    plt.savefig('hubble_resid_variants.png',dpi=200)
    import pdb; pdb.set_trace()
    
if __name__ == "__main__":
    main()
