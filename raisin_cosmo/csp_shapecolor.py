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
_goodcids = np.append(_goodcids,np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',
                                           unpack=True,dtype=str))
_goodcids = np.append(_goodcids,np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',
                                           unpack=True,dtype=str))

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

def runall_fullset():
    #

    dlmags,dlmagerrs,hostmasses,snids,zs,sbvs,sbverrs,avs,dlmags_nocorr,dlmagerrs_nocorr = \
        np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
    
    for survey in ['PS1','DES','CSP']:
        fro = txtobj(f'output/fit_optical/{survey}_RAISIN_optnir.FITRES.TEXT',fitresheader=True)

        for j,i in enumerate(fro.CID):
            if i not in _goodcids: continue

            av_single = fro.AV[j]
            st_single = fro.STRETCH[j]
            sbverrs = np.append(sbverrs,fro.STRETCHERR[j])
            
            with open(f'fit/{survey}_RAISIN_shapecolor.nml') as fin, \
                 open(f'fit/{survey}_RAISIN_shapecolor_tmp.nml','w') as fout:

                for line in fin:
                    #if 'INIVAL_AV' in line:
                    #    print(f'         INIVAL_AV = {av_single:.3f}',file=fout)
                    if 'INIVAL_SHAPE' in line:
                        print(f'         INIVAL_SHAPE = {st_single:.3f}',file=fout)
                    elif 'TEXTFILE_PREFIX' in line:
                        print(f"         TEXTFILE_PREFIX  = 'output/fit_nir/{survey}_RAISIN_NIR_SHAPECOLOR_TMP'",file=fout)
                    elif 'SNCCID_LIST' in line:
                        print(f"         SNCCID_LIST = '{i}'",file=fout)
                    elif 'INISTP_AV' in line:
                        print(f"         INISTP_AV = 0",file=fout)
                    elif 'INISTP_SHAPE' in line:
                        print(f"         INISTP_SHAPE = 0",file=fout)

                    else:
                        print(line.replace('\n',''),file=fout)

            os.system(f'snlc_fit.exe fit/{survey}_RAISIN_shapecolor_tmp.nml')
            frt = txtobj(f'output/fit_nir/{survey}_RAISIN_NIR_SHAPECOLOR_TMP.FITRES.TEXT',fitresheader=True)
            dlmags = np.append(dlmags,frt.DLMAG[0])
            dlmagerrs = np.append(dlmagerrs,frt.DLMAGERR[0])
            zs = np.append(zs,frt.zHD[0])
            snids = np.append(snids,frt.CID[0])
            hostmasses = np.append(hostmasses,frt.HOST_LOGMASS[0])
            sbvs = np.append(sbvs,frt.STRETCH[0])
            avs = np.append(avs,frt.AV[0])

            os.system(f'snlc_fit.exe fit/{survey}_RAISIN_shapecolor_tmp.nml INIVAL_SHAPE 1.0')
            frt = txtobj(f'output/fit_nir/{survey}_RAISIN_NIR_SHAPECOLOR_TMP.FITRES.TEXT',fitresheader=True)
            dlmags_nocorr = np.append(dlmags_nocorr,frt.DLMAG[0])
            dlmagerrs_nocorr = np.append(dlmagerrs_nocorr,frt.DLMAGERR[0])
            #import pdb; pdb.set_trace()
            
    print(np.std(dlmags-cosmo.mu(zs)))
    with open('st_av_corr_mags.txt','w') as fout:
        print('# CID resid DLMAGERR zHD HOST_LOGMASS STRETCH STRETCHERR AV resid2 DLMAGERR2',file=fout)
        for d,de,z,s,h,sbv,se,av,d2,de2 in zip(dlmags,dlmagerrs,zs,snids,hostmasses,sbvs,sbverrs,avs,dlmags_nocorr,dlmagerrs_nocorr):
            print(f"{s} {d-cosmo.mu(z)} {de} {z} {h} {sbv} {se} {av} {d2-cosmo.mu(z)} {de2}",file=fout)

    #0.1421244699193456 - shape from NIR and AV from opt
    #0.16124328730844295 - shape = 1 and AV from opt
    #default: 0.210
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
    #optvnir_stretch()
    runall_fullset()
