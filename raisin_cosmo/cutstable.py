#!/usr/bin/env python
import numpy as np
import os
import snana
from txtobj import txtobj
import cosmo
from astropy.stats import sigma_clipped_stats

# $SNANA_DIR/util/plot_snana.py -v CSPDR3_RAISIN -a -20 -f $RAISIN_ROOT/cosmo/fit/CSP_RAISIN_optnir.nml --plotAll

_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_CIDS.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/PS1_CIDS.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/DES_CIDS.LIST',unpack=True,dtype=str)))
bad_host_list = ['DES16E2cxw','DES16E2rd','DES16X1cpf','PScA470041','PScK450082']
nonia_list = ['PScF520107','2006bt']
lcp = txtobj('output/fit_optical/CSP_RAISIN_optnir.LCPLOT.TEXT',fitresheader=True,rowprfx='OBS') #'output/cutslists/CSP_RAISIN_LCPLOT.TEXT',fitresheader=True,rowprfx='OBS')

#2007ai
# need to leave out 09le

def apply_all_cuts(fr,fropt,restrict_to_good_list=False,stmin=0.8,stmax=1.3):

    # AV
    iGoodAV = np.zeros(len(fr.CID),dtype=bool)
    for j,i in enumerate(fr.CID):
        if i in fropt.CID and fropt.AV[fropt.CID == i] < 0.3*fropt.RV[fropt.CID == i]:
            iGoodAV[j] = True

    # reasonable stretch
    iGoodSt = np.zeros(len(fr.CID),dtype=bool)
    for j,i in enumerate(fr.CID):
        if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > stmin and fropt.STRETCH[fropt.CID == i] < stmax:
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


def main(write=False):

    avmax,stmin,stmax = 1.0,0.75,1.185 #31 #0.99,0.77,1.31

    cspnirfitres = os.path.expandvars('$RAISIN_ROOT/cosmo/output/fit_nir_sys/CSP_RAISIN/CSPDR3_RAISIN/FITOPT000.FITRES.gz')
    ps1nirfitres = os.path.expandvars('$RAISIN_ROOT/cosmo/output/fit_nir_sys/PS1_RAISIN/PS1_RAISIN/FITOPT000.FITRES.gz')
    desnirfitres = os.path.expandvars('$RAISIN_ROOT/cosmo/output/fit_nir_sys/DES_RAISIN/DES_RAISIN/FITOPT000.FITRES.gz')
    frcspnir = txtobj(desnirfitres,fitresheader=True)
    frps1nir = txtobj(desnirfitres,fitresheader=True)
    frdesnir = txtobj(desnirfitres,fitresheader=True)
    cspresid = frcspnir.DLMAG-cosmo.mu(frcspnir.zHD)
    ps1resid = frps1nir.DLMAG-cosmo.mu(frps1nir.zHD)
    desresid = frdesnir.DLMAG-cosmo.mu(frdesnir.zHD)
    cspresid -= np.median(cspresid); ps1resid -= np.median(ps1resid); desresid -= np.median(desresid)
    resid = np.concatenate((cspresid,ps1resid,desresid))
    allcid = np.concatenate((frcspnir.CID,frps1nir.CID,frdesnir.CID))
    residstd = sigma_clipped_stats(resid)[2]; residmean = np.mean(resid)
    chauvenet_list = np.unique(allcid[np.where(abs(resid) > 2.7*residstd)[0]])
    
    #cuts:
    # data quality
    # classification
    # redshift
    # tmax
    # AV
    # stretch

    # ID sighost class z tmax AV st

    # put together CSP list
    cspfiles = np.loadtxt('data/Photometry/bkp/CSPDR3_RAISIN_210426/CSPDR3_RAISIN_FULL.LIST',unpack=True,dtype=str)
    cspfitres = os.path.expandvars('$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT')
    cspids,cspzs,csptpk,csptpk2,cspnir = np.array([]),np.array([]),np.array([],dtype=bool),np.array([],dtype=bool),np.array([],dtype=bool)
    for c in cspfiles:
        sn = snana.SuperNova(f'data/Photometry/bkp/CSPDR3_RAISIN_210426/{c}')
        cspids = np.append(cspids,sn.SNID)
        cspzs = np.append(cspzs,float(sn.REDSHIFT_FINAL.split('+-')[0]))
        csptpk = np.append(csptpk,np.min(sn.MJD) < float(sn.PEAKMJD.split()[0]))
        if len(lcp.Tobs[lcp.CID == sn.SNID]):
            csptpk2 = np.append(csptpk2,np.min(lcp.Tobs[lcp.CID == sn.SNID]) < 0)
        else:
            csptpk2 = np.append(csptpk2,False)
        #if sn.SNID == '2004gs': import pdb; pdb.set_trace()
        if '2009le' not in c:
            cspnir = np.append(cspnir,len(sn.FLT[(sn.FLT == 'Y') | (sn.FLT == 'y') | (sn.FLT == 'J') |
                                                 (sn.FLT == 'j') | (sn.FLT == 'H') | (sn.FLT == 'h')]) > 0)
        else:
            cspnir = np.append(cspnir,False)
    csp_ntot = len(cspids)
    csp_classcut,csp_classcut_ncut = len(cspids)-1,1
    csp_nzcut,csp_nzcut_ncut = len(cspids[(cspzs > 0.01) & (cspids != '2006bt')]),len(cspids)-len(cspids[(cspzs > 0.01) & (cspids != '2006bt')])
    csp_ntpkcut,csp_ntpkcut_ncut = len(cspids[(cspzs > 0.01) & csptpk & (cspids != '2006bt')]),\
        len(cspids[(cspzs > 0.01) & (cspids != '2006bt')])-len(cspids[(cspzs > 0.01) & csptpk & (cspids != '2006bt')])
    csp_ntpkcut2,csp_ntpkcut_ncut2 = len(cspids[(cspzs > 0.01) & csptpk & csptpk2 & (cspids != '2006bt')]),\
        len(cspids[(cspzs > 0.01) & csptpk & (cspids != '2006bt')])-len(cspids[(cspzs > 0.01) & csptpk & csptpk2 & (cspids != '2006bt')])
    csp_nnircut,csp_nnircut_ncut = len(cspids[(cspzs > 0.01) & csptpk & csptpk2 & cspnir & (cspids != '2006bt')]),\
        len(cspids[(cspzs > 0.01) & csptpk2 & csptpk & (cspids != '2006bt')])-len(cspids[(cspzs > 0.01) & csptpk & csptpk2 & cspnir & (cspids != '2006bt')])
    cspids_prefr_final = cspids[(cspzs > 0.01) & csptpk & csptpk2 & cspnir & (cspids != '2006bt')]

    frcsp = txtobj(cspfitres,fitresheader=True)
    csptpkfr = np.array([],dtype=bool)
    csptpkofr = np.array([],dtype=bool)
    cspnirfr = np.array([],dtype=bool)
    for j,i in enumerate(frcsp.CID):
        sn = snana.SuperNova(f'data/Photometry/bkp/CSPDR3_RAISIN_210426/CSPDR3_{i}.PKMJD.DAT')
        #if len(lcp.Tobs[lcp.CID == i]):
        #    csptpkfr = np.append(csptpkfr,np.min(lcp.Tobs[(lcp.CID == i) & (lcp.DATAFLAG == 1)]) < 0)
        #else:
        #    csptpkfr = np.append(csptpkfr,False)
        csptpkofr = np.append(csptpkofr,np.min(sn.MJD) < float(sn.PEAKMJD.split()[0]))
#        import pdb; pdb.set_trace()
        cspnirfr = np.append(cspnirfr,len(sn.FLT[((sn.FLT == 'Y') | (sn.FLT == 'y') | (sn.FLT == 'J') | (sn.FLT == 'j') | (sn.FLT == 'H') | (sn.FLT == 'h')) & (sn.MJD-float(sn.PEAKMJD.split()[0]) > -5)]) > 0)

    cspids_fr_final2 = frcsp.CID[(frcsp.zCMB > 0.01) & cspnirfr & (frcsp.CID != '2006bt')]
    cspids_fr_final = frcsp.CID[(frcsp.zCMB > 0.01) & cspnirfr & csptpkofr & (frcsp.CID != '2006bt')]
    for i in cspids_fr_final2:
        if i not in cspids_fr_final: print(i)
    #import pdb; pdb.set_trace()
    iCut = ((frcsp.zCMB > 0.01) & csptpkofr & cspnirfr) | (frcsp.CID == '2006bt') #csptpkfr & 
    for k in frcsp.__dict__.keys():
        frcsp.__dict__[k] = frcsp.__dict__[k][iCut]

    csp_fitcut,csp_fitcut_ncut = len(cspids_fr_final),len(cspids_prefr_final)-len(cspids_fr_final)
    csp_navcut,csp_navcut_ncut = len(frcsp.CID[frcsp.AV < avmax]),len(frcsp.CID[frcsp.AV > avmax])
    csp_nstcut,csp_nstcut_ncut = len(frcsp.CID[(frcsp.AV < avmax) & (frcsp.STRETCH > stmin) & (frcsp.STRETCH < stmax)]),\
                                 len(frcsp.CID[frcsp.AV < avmax])-len(frcsp.CID[(frcsp.AV < avmax) & (frcsp.STRETCH > stmin) & (frcsp.STRETCH < stmax)])
    csp_nsterrcut,csp_nsterrcut_ncut = \
        len(frcsp.CID[(frcsp.AV < avmax) & (frcsp.STRETCH > stmin) & (frcsp.STRETCH < stmax) & (frcsp.STRETCHERR < 0.3)]),\
        len(frcsp.CID[(frcsp.AV < avmax) & (frcsp.STRETCH > stmin) & (frcsp.STRETCH < stmax)])-\
        len(frcsp.CID[(frcsp.AV < avmax) & (frcsp.STRETCH > stmin) & (frcsp.STRETCH < stmax) & (frcsp.STRETCHERR < 0.3)])
    import pdb; pdb.set_trace()

    iGood = np.array([],dtype=int)
    iChauv = np.array([],dtype=int)
    for j,i in enumerate(frcsp.CID):
        if i not in chauvenet_list:
            iGood = np.append(iGood,j)
        elif i in chauvenet_list:
            iChauv = np.append(iChauv,j)
    for k in frcsp.__dict__.keys():
        frcsp.__dict__[k] = frcsp.__dict__[k][iGood]

    csp_nchauvenetcut = \
        len(frcsp.CID[(frcsp.AV < avmax) & (frcsp.STRETCH > stmin) & (frcsp.STRETCH < stmax) & (frcsp.STRETCHERR < 0.3)])
    cspnchauvenetcut_ncut = len(iChauv)

    

    ps1files,ps1cids = np.loadtxt('data/Photometry/PS1_RAISIN/PS1_RAISIN_FULL.LIST',unpack=True,dtype=str)
    ps1fitres = os.path.expandvars('$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT')

    ps1_ndqcut_ncut = 0
    ps1_classcut_ncut = 0
    for p in ps1cids:
        if p in bad_host_list:
            ps1_ndqcut_ncut += 1
        elif p in nonia_list:
            ps1_classcut_ncut += 1
    
    ps1_ndqcut = len(ps1files)-ps1_ndqcut_ncut
    ps1_classcut = len(ps1files)-ps1_ndqcut_ncut-ps1_classcut_ncut
    frps1 = txtobj(ps1fitres,fitresheader=True)
    iGood = np.array([],dtype=int)
    for j,i in enumerate(frps1.CID):
        if i not in bad_host_list and i not in nonia_list:
            iGood = np.append(iGood,j)
    for k in frps1.__dict__.keys():
        frps1.__dict__[k] = frps1.__dict__[k][iGood]
    
    ps1_navcut,ps1_navcut_ncut = len(frps1.CID[frps1.AV < avmax]),len(frps1.CID[frps1.AV > avmax])
    ps1_nstcut,ps1_nstcut_ncut = len(frps1.CID[(frps1.AV < avmax) & (frps1.STRETCH > stmin) & (frps1.STRETCH < stmax)]),\
                                 len(frps1.CID[frps1.AV < avmax])-len(frps1.CID[(frps1.AV < avmax) & (frps1.STRETCH > stmin) & (frps1.STRETCH < stmax)])
    ps1_nsterrcut,ps1_nsterrcut_ncut = \
        len(frps1.CID[(frps1.AV < avmax) & (frps1.STRETCH > stmin) & (frps1.STRETCH < stmax) & (frps1.STRETCHERR < 0.3)]),\
        len(frps1.CID[(frps1.AV < avmax) & (frps1.STRETCH > stmin) & (frps1.STRETCH < stmax)])-\
        len(frps1.CID[(frps1.AV < avmax) & (frps1.STRETCH > stmin) & (frps1.STRETCH < stmax) & (frps1.STRETCHERR < 0.3)])

    iGood = np.array([],dtype=int)
    iChauv = np.array([],dtype=int)
    for j,i in enumerate(frps1.CID):
        if i not in chauvenet_list:
            iGood = np.append(iGood,j)
        elif i in chauvenet_list:
            iChauv = np.append(iChauv,j)
    for k in frps1.__dict__.keys():
        frps1.__dict__[k] = frps1.__dict__[k][iGood]

    ps1_nchauvenetcut = \
        len(frps1.CID[(frps1.AV < avmax) & (frps1.STRETCH > stmin) & (frps1.STRETCH < stmax) & (frps1.STRETCHERR < 0.3)])
    ps1nchauvenetcut_ncut = len(iChauv)

    
    desfiles,descids = np.loadtxt('data/Photometry/DES_RAISIN/DES_RAISIN_FULL.LIST',unpack=True,dtype=str)
    desfitres = os.path.expandvars('$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT')
    
    des_ndqcut_ncut = 0
    des_classcut_ncut = 0
    for d in descids:
        if d in bad_host_list:
            des_ndqcut_ncut += 1
        elif d in nonia_list:
            des_classcut_ncut += 1
    
    des_ndqcut = len(desfiles)-des_ndqcut_ncut
    des_classcut = len(desfiles)-des_ndqcut_ncut-des_classcut_ncut
    frdes = txtobj(desfitres,fitresheader=True)
    iGood = np.array([],dtype=int)
    for j,i in enumerate(frdes.CID):
        if i not in bad_host_list and i not in nonia_list:
            iGood = np.append(iGood,j)
    for k in frdes.__dict__.keys():
        frdes.__dict__[k] = frdes.__dict__[k][iGood]


    des_navcut,des_navcut_ncut = len(frdes.CID[frdes.AV < avmax]),len(frdes.CID[frdes.AV > avmax])
    des_nstcut,des_nstcut_ncut = len(frdes.CID[(frdes.AV < avmax) & (frdes.STRETCH > stmin) & (frdes.STRETCH < stmax)]),\
                                 len(frdes.CID[frdes.AV < avmax])-len(frdes.CID[(frdes.AV < avmax) & (frdes.STRETCH > stmin) & (frdes.STRETCH < stmax)])
    des_nsterrcut,des_nsterrcut_ncut = \
        len(frdes.CID[(frdes.AV < avmax) & (frdes.STRETCH > stmin) & (frdes.STRETCH < stmax) & (frdes.STRETCHERR < 0.3)]),\
        len(frdes.CID[(frdes.AV < avmax) & (frdes.STRETCH > stmin) & (frdes.STRETCH < stmax)])-\
        len(frdes.CID[(frdes.AV < avmax) & (frdes.STRETCH > stmin) & (frdes.STRETCH < stmax) & (frdes.STRETCHERR < 0.3)])



    iGood = np.array([],dtype=int)
    iChauv = np.array([],dtype=int)
    for j,i in enumerate(frdes.CID):
        if i not in chauvenet_list:
            iGood = np.append(iGood,j)
        elif i in chauvenet_list:
            iChauv = np.append(iChauv,j)
    for k in frdes.__dict__.keys():
        frdes.__dict__[k] = frdes.__dict__[k][iGood]

    des_nchauvenetcut = \
        len(frdes.CID[(frdes.AV < avmax) & (frdes.STRETCH > stmin) & (frdes.STRETCH < stmax) & (frdes.STRETCHERR < 0.3)])
    desnchauvenetcut_ncut = len(iChauv)
    
    
    initline = f"&\\nodata&{csp_ntot:.0f}&\\nodata&{len(ps1files):.0f}&\\nodata&{len(desfiles):.0f}\\\\"
    dqline = f"$\\sigma_{{host}} < 0.1$ mag&0&{csp_ntot}&{ps1_ndqcut_ncut}&{ps1_ndqcut:.0f}&{des_ndqcut_ncut}&{des_ndqcut:.0f}\\\\"
    classline = f"normal SN\\,Ia&{csp_classcut_ncut}&{csp_classcut:.0f}&{ps1_classcut_ncut}&{ps1_classcut:.0f}&{des_classcut_ncut}&{des_classcut:.0f}\\\\"
    zline = f"$z > 0.01$&{csp_nzcut_ncut:.0f}&{csp_nzcut:.0f}&0&{ps1_classcut:.0f}&0&{des_classcut:.0f}\\\\"
    tmaxline = f"pre-max&{csp_ntpkcut_ncut:.0f}&{csp_ntpkcut:.0f}&0&{ps1_classcut:.0f}&0&{des_classcut:.0f}\\\\"
    tmaxline2 = f"pre-max, $\\chi^2$ cut&{csp_ntpkcut_ncut2:.0f}&{csp_ntpkcut2:.0f}&0&{ps1_classcut:.0f}&0&{des_classcut:.0f}\\\\"
    nirline = f"NIR data&{csp_nnircut_ncut:.0f}&{csp_nnircut:.0f}&0&{ps1_classcut:.0f}&0&{des_classcut:.0f}\\\\"
    fitline = f"SNooPy fitting&{csp_fitcut_ncut:.0f}&{csp_fitcut:.0f}&0&{ps1_classcut:.0f}&0&{des_classcut:.0f}\\\\"
    avline = f"$A_V < 1.0$&{csp_navcut_ncut}&{csp_navcut}&{ps1_navcut_ncut}&{ps1_navcut:.0f}&{des_navcut_ncut}&{des_navcut:.0f}\\\\"
    stline = f"$0.75 < s_{{BV}} < 1.18$&{csp_nstcut_ncut}&{csp_nstcut}&{ps1_nstcut_ncut}&{ps1_nstcut:.0f}&{des_nstcut_ncut}&{des_nstcut:.0f}\\\\"
    sterrline = f"$\sigma(s_{{BV}}) < 0.2$&{csp_nsterrcut_ncut}&{csp_nsterrcut}&{ps1_nsterrcut_ncut}&{ps1_nsterrcut:.0f}&{des_nsterrcut_ncut}&{des_nsterrcut:.0f}\\\\"
    chauvenetline = f"Chauvenet&{cspnchauvenetcut_ncut}&{csp_nchauvenetcut}&{ps1nchauvenetcut_ncut}&{ps1_nchauvenetcut}&{desnchauvenetcut_ncut}&{des_nchauvenetcut}\\\\"

    #import pdb; pdb.set_trace()
    for line in [initline,dqline,classline,zline,tmaxline,tmaxline2,nirline,fitline,avline,stline,sterrline,chauvenetline]:
        print(line)

    if write:
        with open('output/goodcids/DES_GOODCIDS_LATEST.LIST','w') as fout:
            for i in frdes.CID[(frdes.AV < avmax) & (frdes.STRETCH > stmin) & (frdes.STRETCH < stmax) & (frdes.STRETCHERR < 0.2)]:
                print(i,file=fout)
        with open('output/goodcids/PS1_GOODCIDS_LATEST.LIST','w') as fout:
            for i in frps1.CID[(frps1.AV < avmax) & (frps1.STRETCH > stmin) & (frps1.STRETCH < stmax) & (frps1.STRETCHERR < 0.2)]:
                print(i,file=fout)
        with open('output/goodcids/CSP_GOODCIDS_LATEST.LIST','w') as fout:
            for i in frcsp.CID[(frcsp.AV < avmax) & (frcsp.STRETCH > stmin) & (frcsp.STRETCH < stmax) & (frcsp.STRETCHERR < 0.2)]:
                print(i,file=fout)

#def write_new_pkmjd():
#    files = glob
            
if __name__ == "__main__":
    main()
