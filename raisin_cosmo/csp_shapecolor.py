#!/usr/bin/env python
# D. Jones - 5/12/21

import numpy as np
from txtobj import txtobj
import pylab as plt
plt.ion()
import cosmo
import os

_goodcids = np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',
                       unpack=True,dtype=str)

def main():
    frnc = txtobj('output/fit_nir_sys/CSP_RAISIN/CSPDR3_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    frc = txtobj('output/fit_nir/CSP_RAISIN_NIR_SHAPE.FITRES.TEXT',
                 fitresheader=True)

    idxnc,idxc = np.array([],dtype=int),np.array([],dtype=int)
    for g in _goodcids:
        #if g == '2009Y': continue
        idxnc = np.append(idxnc,np.where(frnc.CID == g)[0])
        idxc = np.append(idxc,np.where(frc.CID == g)[0])

    for k in frnc.__dict__.keys():
        frnc.__dict__[k] = frnc.__dict__[k][idxnc]
    for k in frc.__dict__.keys():
        frc.__dict__[k] = frc.__dict__[k][idxc]

    for i in frnc.CID:
        if i not in frc.CID: print(i)

    residc = frc.DLMAG - cosmo.mu(frc.zHD)
    residnc = frnc.DLMAG - cosmo.mu(frnc.zHD)
    print(np.std(residc),np.std(residnc))
    # 0.14215787140082425 0.15577450111219995
    # 0.14311934697087422 0.15703609248383055
    import pdb; pdb.set_trace()

def avrun():
    #
    
    fro = txtobj('output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',fitresheader=True)

    dlmags,zs = np.array([]),np.array([])
    for j,i in enumerate(fro.CID):
        if i not in _goodcids: continue

        av_single = fro.AV[j]

        with open('fit/CSP_RAISIN_shapecolor.nml') as fin, \
             open('fit/CSP_RAISIN_shapecolor_tmp.nml','w') as fout:

            for line in fin:
                if 'INIVAL_AV' in line:
                    print(f'         INIVAL_AV = {av_single:.3f}',file=fout)
                elif 'TEXTFILE_PREFIX' in line:
                    print(f"         TEXTFILE_PREFIX  = 'output/fit_nir/CSP_RAISIN_NIR_SHAPECOLOR_TMP'",file=fout)
                elif 'SNCCID_LIST' in line:
                    print(f"         SNCCID_LIST = '{i}'",file=fout)
                elif 'INISTP_AV' in line:
                    print(f"         INISTP_AV = 0",file=fout)
                elif 'INISTP_SHAPE' in line:
                    print(f"         INISTP_SHAPE = 0",file=fout)

                else:
                    print(line.replace('\n',''),file=fout)

        os.system('snlc_fit.exe fit/CSP_RAISIN_shapecolor_tmp.nml')
        frt = txtobj('output/fit_nir/CSP_RAISIN_NIR_SHAPECOLOR_TMP.FITRES.TEXT',fitresheader=True)
        dlmags = np.append(dlmags,frt.DLMAG[0])
        zs = np.append(zs,frt.zHD[0])

    print(np.std(dlmags-cosmo.mu(zs)))
    #0.1421244699193456 - shape from NIR and AV from opt
    #0.16124328730844295 - shape = 1 and AV from opt
    import pdb; pdb.set_trace()

def optvnir_stretch():

    fro = txtobj('output/fit_optical/CSP_RAISIN_optical.FITRES.TEXT',fitresheader=True)
    frn = txtobj('output/fit_nir/CSP_RAISIN_NIR_SHAPE.FITRES.TEXT',fitresheader=True)

    sbv_opt,sbv_opterr,sbv_nir,sbv_nirerr = np.array([]),np.array([]),np.array([]),np.array([])
    for j,i in enumerate(fro.CID):
        if i not in _goodcids: continue
        if i not in frn.CID: continue

        if fro.STRETCHERR[j] > 0.2: continue
        if frn.STRETCHERR[frn.CID == i][0] > 0.2: continue
        if fro.STRETCH[j] > 1.3: continue
        if frn.STRETCH[frn.CID == i][0] > 1.3: continue
        sbv_opt = np.append(sbv_opt,fro.STRETCH[j])
        sbv_nir = np.append(sbv_nir,frn.STRETCH[frn.CID == i][0])
        sbv_opterr = np.append(sbv_opterr,fro.STRETCHERR[j])
        sbv_nirerr = np.append(sbv_nirerr,frn.STRETCHERR[frn.CID == i][0])

    sbv = np.linspace(0.8,1.3,100)
    plt.plot(sbv,sbv,color='k')
    plt.errorbar(sbv_opt,sbv_nir,xerr=sbv_opterr,yerr=sbv_nirerr,fmt='.')
    plt.xlabel('$s_{BV}$ from optical data')
    plt.ylabel('$s_{BV}$ from NIR data')
    plt.xlim([0.8,1.22])
    plt.ylim([0.8,1.22])
    import pdb; pdb.set_trace()
    
if __name__ == "__main__":
    #main()
    #avrun()
    optvnir_stretch()
