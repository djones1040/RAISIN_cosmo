#!/usr/bin/env python
# D. Jones 2/15/21

import numpy as np
import snana
import glob
from coordutils import GetSexigesimalString

def main():

    
    for files in [glob.glob('data/photometry/PS1_RAISIN/*.snana.dat'),
                  glob.glob('data/photometry/DES_RAISIN/*.snana.dat')]:
        zlist = []
        for f in files:
            sn = snana.SuperNova(f)
            zlist += [float(sn.REDSHIFT_FINAL.split()[0])]

        idx = np.argsort(zlist)
        for f in np.array(files)[idx]:
            sn = snana.SuperNova(f)
            if 'DECL' in sn.__dict__.keys():
                ra_str,dec_str = GetSexigesimalString(sn.RA,sn.DECL)
            else:
                ra_str,dec_str = GetSexigesimalString(sn.RA,sn.DEC)
            first_mjd = np.min(sn.MJD[sn.FLUXCAL/sn.FLUXCALERR > 5])
            print(f"{sn.SNID}&{ra_str}&{dec_str}&{float(sn.REDSHIFT_FINAL.split()[0]):.3f}&{first_mjd:.0f}\\\\")
        
if __name__ == "__main__":
    main()
