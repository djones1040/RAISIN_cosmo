#!/usr/bin/env python
import numpy as np
import os
import snana
from txtobj import txtobj

# $SNANA_DIR/util/plot_snana.py -v CSPDR3_RAISIN -a -20 -f $RAISIN_ROOT/cosmo/fit/CSP_RAISIN_optnir.nml --plotAll

_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_CIDS.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/PS1_CIDS.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/DES_CIDS.LIST',unpack=True,dtype=str)))
bad_host_list = ['DES16E2cxw','DES16E2rd','DES16X1cpf','PScA470041','PScK450082']
nonia_list = ['PScF520107']
lcp = txtobj('output/cutslists/CSP_RAISIN_LCPLOT.TEXT',fitresheader=True,rowprfx='OBS')

def apply_all_cuts(fr,fropt,restrict_to_good_list=False):

    # AV
    iGoodAV = np.zeros(len(fr.CID),dtype=bool)
    for j,i in enumerate(fr.CID):
        if i in fropt.CID and fropt.AV[fropt.CID == i] < 0.3*fropt.RV[fropt.CID == i]:
            iGoodAV[j] = True

    # reasonable stretch
    iGoodSt = np.zeros(len(fr.CID),dtype=bool)
    for j,i in enumerate(fr.CID):
        if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > 0.8 and fropt.STRETCH[fropt.CID == i] < 1.3:
            iGoodSt[j] = True

    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iGoodAV & iGoodSt]

    if restrict_to_good_list:
        iGood = np.array([],dtype=int)
        for j,i in enumerate(fr.CID):
            if i in _goodcids: iGood = np.append(iGood,j)
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]

    return fr


def main():

    #cuts:
    # data quality
    # classification
    # redshift
    # tmax
    # AV
    # stretch

    # ID sighost class z tmax AV st

    # put together CSP list
    cspfiles = np.loadtxt('data/Photometry/CSPDR3_RAISIN/CSPDR3_RAISIN.LIST',unpack=True,dtype=str)
    cspfitres = os.path.expandvars('$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT')
    cspids,cspzs,csptpk,csptpk2,cspnir = np.array([]),np.array([]),np.array([],dtype=bool),np.array([],dtype=bool),np.array([],dtype=bool)
    for c in cspfiles:
        sn = snana.SuperNova(f'data/Photometry/CSPDR3_RAISIN/{c}')
        cspids = np.append(cspids,sn.SNID)
        cspzs = np.append(cspzs,float(sn.REDSHIFT_FINAL.split('+-')[0]))
        csptpk = np.append(csptpk,np.min(sn.MJD) < float(sn.PEAKMJD.split()[0]))
        if len(lcp.Tobs[lcp.CID == sn.SNID]):
            csptpk2 = np.append(csptpk2,np.min(lcp.Tobs[lcp.CID == sn.SNID]) < 0)
        else:
            csptpk2 = np.append(csptpk2,False)
        cspnir = np.append(cspnir,len(sn.FLT[(sn.FLT == 'Y') | (sn.FLT == 'y') | (sn.FLT == 'J') | (sn.FLT == 'j') | (sn.FLT == 'H') | (sn.FLT == 'h')]) > 0)
    csp_ntot = len(cspids)
    csp_nzcut,csp_nzcut_ncut = len(cspids[cspzs > 0.01]),len(cspids)-len(cspids[cspzs > 0.01])
    csp_ntpkcut,csp_ntpkcut_ncut = len(cspids[(cspzs > 0.01) & csptpk]),len(cspids[cspzs > 0.01])-len(cspids[(cspzs > 0.01) & csptpk])
    csp_ntpkcut2,csp_ntpkcut_ncut2 = len(cspids[(cspzs > 0.01) & csptpk & csptpk2]),len(cspids[(cspzs > 0.01) & csptpk])-len(cspids[(cspzs > 0.01) & csptpk & csptpk2])
    csp_nnircut,csp_nnircut_ncut = len(cspids[(cspzs > 0.01) & csptpk & csptpk2 & cspnir]),len(cspids[(cspzs > 0.01) & csptpk2 & csptpk])-len(cspids[(cspzs > 0.01) & csptpk & csptpk2 & cspnir])
    cspids_prefr_final = cspids[(cspzs > 0.01) & csptpk & csptpk2 & cspnir]
                                
    frcsp = txtobj(cspfitres,fitresheader=True)
    csptpkfr = np.array([],dtype=bool)
    csptpkofr = np.array([],dtype=bool)
    cspnirfr = np.array([],dtype=bool)
    for j,i in enumerate(frcsp.CID):
        sn = snana.SuperNova(f'data/Photometry/CSPDR3_RAISIN/CSPDR3_{i}.PKMJD.DAT')
        if len(lcp.Tobs[lcp.CID == i]):
            csptpkfr = np.append(csptpkfr,np.min(lcp.Tobs[(lcp.CID == i) & (lcp.DATAFLAG == 1)]) < 0)
        else:
            csptpkfr = np.append(csptpkfr,False)
        csptpkofr = np.append(csptpkofr,np.min(sn.MJD) < float(sn.PEAKMJD.split()[0]))
        cspnirfr = np.append(cspnirfr,len(sn.FLT[(sn.FLT == 'Y') | (sn.FLT == 'y') | (sn.FLT == 'J') | (sn.FLT == 'j') | (sn.FLT == 'H') | (sn.FLT == 'h')]) > 0)

    cspids_fr_final2 = frcsp.CID[csptpkfr & (frcsp.zCMB > 0.01) & cspnirfr & csptpkfr]
    cspids_fr_final = frcsp.CID[csptpkfr & (frcsp.zCMB > 0.01) & cspnirfr & csptpkofr & csptpkfr]
    #for i in cspids_prefr_final:
    #    if i not in cspids_fr_final: print(i)
    for i in cspids_fr_final2:
        if i not in cspids_fr_final: print(i)
    import pdb; pdb.set_trace()
    iCut = csptpkfr & (frcsp.zCMB > 0.01) & csptpkofr & cspnirfr
    for k in frcsp.__dict__.keys():
        frcsp.__dict__[k] = frcsp.__dict__[k][iCut]

    csp_fitcut,csp_fitcut_ncut = len(cspids_fr_final),len(cspids_prefr_final)-len(cspids_fr_final)
    csp_navcut,csp_navcut_ncut = len(frcsp.CID[frcsp.AV < 1.0]),len(frcsp.CID[frcsp.AV > 1.0])
    csp_nstcut,csp_nstcut_ncut = len(frcsp.CID[(frcsp.AV < 1.0) & (frcsp.STRETCH > 0.8) & (frcsp.STRETCH < 1.3)]),\
                                 len(frcsp.CID[frcsp.AV < 1.0])-len(frcsp.CID[(frcsp.AV < 1.0) & (frcsp.STRETCH > 0.8) & (frcsp.STRETCH < 1.3)])
    csp_nsterrcut,csp_nsterrcut_ncut = \
        len(frcsp.CID[(frcsp.AV < 1.0) & (frcsp.STRETCH > 0.8) & (frcsp.STRETCH < 1.3) & (frcsp.STRETCHERR < 0.3)]),\
        len(frcsp.CID[(frcsp.AV < 1.0) & (frcsp.STRETCH > 0.8) & (frcsp.STRETCH < 1.3)])-\
        len(frcsp.CID[(frcsp.AV < 1.0) & (frcsp.STRETCH > 0.8) & (frcsp.STRETCH < 1.3) & (frcsp.STRETCHERR < 0.3)])

    ps1files = np.loadtxt('data/Photometry/PS1_RAISIN/PS1_RAISIN_FULL.LIST',unpack=True,dtype=str)
    ps1fitres = os.path.expandvars('$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT')
    ps1_ndqcut,ps1_ndqcut_ncut = len(ps1files)-2,2
    ps1_classcut,ps1_classcut_ncut = len(ps1files)-3,1
    frps1 = txtobj(ps1fitres,fitresheader=True)
    iGood = frps1.CID != 'PScA470041'
    for k in frps1.__dict__.keys():
        frps1.__dict__[k] = frps1.__dict__[k][iGood]
    
    ps1_navcut,ps1_navcut_ncut = len(frps1.CID[frps1.AV < 1.0]),len(frps1.CID[frps1.AV > 1.0])
    ps1_nstcut,ps1_nstcut_ncut = len(frps1.CID[(frps1.AV < 1.0) & (frps1.STRETCH > 0.8) & (frps1.STRETCH < 1.3)]),\
                                 len(frps1.CID[frps1.AV < 1.0])-len(frps1.CID[(frps1.AV < 1.0) & (frps1.STRETCH > 0.8) & (frps1.STRETCH < 1.3)])
    ps1_nsterrcut,ps1_nsterrcut_ncut = \
        len(frps1.CID[(frps1.AV < 1.0) & (frps1.STRETCH > 0.8) & (frps1.STRETCH < 1.3) & (frps1.STRETCHERR < 0.3)]),\
        len(frps1.CID[(frps1.AV < 1.0) & (frps1.STRETCH > 0.8) & (frps1.STRETCH < 1.3)])-\
        len(frps1.CID[(frps1.AV < 1.0) & (frps1.STRETCH > 0.8) & (frps1.STRETCH < 1.3) & (frps1.STRETCHERR < 0.3)])
    
    desfiles = np.loadtxt('data/Photometry/DES_RAISIN/DES_RAISIN_FULL.LIST',unpack=True,dtype=str)
    desfitres = os.path.expandvars('$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT')
    des_ndqcut,des_ndqcut_ncut = len(desfiles)-3,3
    des_classcut,des_classcut_ncut = len(desfiles)-3,0
    frdes = txtobj(desfitres,fitresheader=True)
    des_navcut,des_navcut_ncut = len(frdes.CID[frdes.AV < 1.0]),len(frdes.CID[frdes.AV > 1.0])
    des_nstcut,des_nstcut_ncut = len(frdes.CID[(frdes.AV < 1.0) & (frdes.STRETCH > 0.8) & (frdes.STRETCH < 1.3)]),\
                                 len(frdes.CID[frdes.AV < 1.0])-len(frdes.CID[(frdes.AV < 1.0) & (frdes.STRETCH > 0.8) & (frdes.STRETCH < 1.3)])
    des_nsterrcut,des_nsterrcut_ncut = \
        len(frdes.CID[(frdes.AV < 1.0) & (frdes.STRETCH > 0.8) & (frdes.STRETCH < 1.3) & (frdes.STRETCHERR < 0.3)]),\
        len(frdes.CID[(frdes.AV < 1.0) & (frdes.STRETCH > 0.8) & (frdes.STRETCH < 1.3)])-\
        len(frdes.CID[(frdes.AV < 1.0) & (frdes.STRETCH > 0.8) & (frdes.STRETCH < 1.3) & (frdes.STRETCHERR < 0.3)])
    
    
    initline = f"&\\nodata&{csp_ntot:.0f}&\\nodata&{len(ps1files):.0f}&\\nodata&{len(desfiles):.0f}\\\\"
    dqline = f"$\\sigma_{{host}} < 0.1$ mag&0&{csp_ntot}&{ps1_ndqcut_ncut}&{ps1_ndqcut:.0f}&{des_ndqcut_ncut}&{des_ndqcut:.0f}\\\\"
    classline = f"confirmed SN\\,Ia&0&{csp_ntot:.0f}&{ps1_classcut_ncut}&{ps1_classcut:.0f}&{des_classcut_ncut}&{des_classcut:.0f}\\\\"
    zline = f"$z > 0.01$&{csp_nzcut_ncut:.0f}&{csp_nzcut:.0f}&0&{ps1_classcut:.0f}&0&{des_classcut:.0f}\\\\"
    tmaxline = f"pre-max&{csp_ntpkcut_ncut:.0f}&{csp_ntpkcut:.0f}&0&{ps1_classcut:.0f}&0&{des_classcut:.0f}\\\\"
    tmaxline2 = f"pre-max, $\\chi^2$ cut&{csp_ntpkcut_ncut2:.0f}&{csp_ntpkcut2:.0f}&0&{ps1_classcut:.0f}&0&{des_classcut:.0f}\\\\"
    nirline = f"NIR data&{csp_nnircut_ncut:.0f}&{csp_nnircut:.0f}&0&{ps1_classcut:.0f}&0&{des_classcut:.0f}\\\\"
    fitline = f"SNooPy fitting&{csp_fitcut_ncut:.0f}&{csp_fitcut:.0f}&0&{ps1_classcut:.0f}&0&{des_classcut:.0f}\\\\"
    avline = f"$A_V < 1.0$&{csp_navcut_ncut}&{csp_navcut}&{ps1_navcut_ncut}&{ps1_navcut:.0f}&{des_navcut_ncut}&{des_navcut:.0f}\\\\"
    stline = f"$0.8 < s_{{BV}} < 1.3$&{csp_nstcut_ncut}&{csp_nstcut}&{ps1_nstcut_ncut}&{ps1_nstcut:.0f}&{des_nstcut_ncut}&{des_nstcut:.0f}\\\\"
    sterrline = f"$\sigma(s_{{BV}}) < 0.2$&{csp_nsterrcut_ncut}&{csp_nsterrcut}&{ps1_nsterrcut_ncut}&{ps1_nsterrcut:.0f}&{des_nsterrcut_ncut}&{des_nsterrcut:.0f}\\\\"
    import pdb; pdb.set_trace()
    for line in [initline,dqline,classline,zline,tmaxline,tmaxline2,nirline,fitline,avline,stline,sterrline]:
        print(line)

    with open('output/goodcids/DES_GOODCIDS_LATEST.LIST','w') as fout:
        for i in frdes.CID[(frdes.AV < 1.0) & (frdes.STRETCH > 0.8) & (frdes.STRETCH < 1.3) & (frdes.STRETCHERR < 0.2)]:
            print(i,file=fout)
    with open('output/goodcids/PS1_GOODCIDS_LATEST.LIST','w') as fout:
        for i in frps1.CID[(frps1.AV < 1.0) & (frps1.STRETCH > 0.8) & (frps1.STRETCH < 1.3) & (frps1.STRETCHERR < 0.2)]:
            print(i,file=fout)
    with open('output/goodcids/CSP_GOODCIDS_LATEST.LIST','w') as fout:
        for i in frcsp.CID[(frcsp.AV < 1.0) & (frcsp.STRETCH > 0.8) & (frcsp.STRETCH < 1.3) & (frcsp.STRETCHERR < 0.2)]:
            print(i,file=fout)

#def write_new_pkmjd():
#    files = glob
            
if __name__ == "__main__":
    main()
