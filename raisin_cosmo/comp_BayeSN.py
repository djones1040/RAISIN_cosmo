#!/usr/bin/env python
# D. Jones - 4/23/21
# see what the typical difference between BayeSN and SNooPy
# distances is for the purpose of systematics

import numpy as np
from txtobj import txtobj

def main():
    frs_lowz = txtobj('output/fit_nir_sys/CSP_RAISIN/CSPDR3_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    frs_raisin = txtobj('output/fit_nir_sys/PS1_RAISIN/PS1_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    frs_des = txtobj('output/fit_nir_sys/DES_RAISIN/DES_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    for k in frs_raisin.__dict__.keys():
        frs_raisin.__dict__[k] = np.append(frs_raisin.__dict__[k],frs_des.__dict__[k])
    
    frb_lowz = txtobj('output/BayeSN/CSP_RAISIN_BayeSN_NIR_theta_-1_dists_jones_snoopy_tmax.txt',fitresheader=False)
    frb_raisin = txtobj('output/BayeSN/RAISIN_BayeSN_theta_-1_NIR_dists_jones_tmax.txt',fitresheader=False)

    # low-z

    muls,mulb = np.array([]),np.array([])
    for j,i in enumerate(frs_lowz.CID):
        iB = np.where(frb_lowz.CID == i)[0]
        muls = np.append(muls,frs_lowz.DLMAG-cosmo.mu(frs_lowz.zHD))
        mulv = np.append(muls,frb_lowz.__dict__['mu_mu+eta'][iB])        
    # raisin
    
    import pdb; pdb.set_trace()
    
if __name__ == "__main__":
    main()
