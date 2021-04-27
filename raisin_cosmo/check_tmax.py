#!/usr/bin/env python
# D. Jones - 4/26/21

import pylab as plt
import numpy as np
from txtobj import txtobj
import snana
import glob

def main():

    frlowz = txtobj('output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',fitresheader=True)
    frps1 = txtobj('output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT',fitresheader=True)
    frdes = txtobj('output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT',fitresheader=True)

    lowzfiles = glob.glob('data/Photometry/CSPDR3_RAISIN/*PKMJD.DAT')
    ps1files = glob.glob('data/Photometry/PS1_RAISIN/*snana.dat')
    desfiles = glob.glob('data/Photometry/DES_RAISIN/*snana.dat')

    lowzcids,ps1cids,descids = \
        'output/goodcids/CSP_GOODCIDS_LATEST.LIST',\
        'output/goodcids/PS1_GOODCIDS_LATEST.LIST',\
        'output/goodcids/DES_GOODCIDS_LATEST.LIST'

    pkmjd_lcfiles,pkmjd_fits,snids = np.array([]),np.array([]),np.array([])
    for files,fr,cids,name in zip(
            [lowzfiles,ps1files,desfiles],
            [frlowz,frps1,frdes],
            [lowzcids,ps1cids,descids],
            ['lowz','PS1','DES']):
        #if name != 'DES': continue
        cidlist = np.loadtxt(cids,unpack=True,dtype=str)
        
        for f in files:
            sn = snana.SuperNova(f)
            if sn.SNID not in cidlist: continue

            try: pkmjd_lcfiles = np.append(pkmjd_lcfiles,float(sn.PEAKMJD.split()[0]))
            except: pkmjd_lcfiles = np.append(pkmjd_lcfiles,sn.PEAKMJD)
            pkmjd_fits = np.append(pkmjd_fits,fr.PKMJD[fr.CID == sn.SNID][0])
            snids = np.append(snids,sn.SNID)
            
    import pdb; pdb.set_trace()

# check 2005hj, 2007as, 2005al
# maybe: 2005al, 2008gl
    
if __name__ == "__main__":
    main()
