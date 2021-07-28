#!/usr/bin/env python

import numpy as np
import os
from txtobj import txtobj

def wfit():
	# run wfit on pantheon with planck+18 data
	# diagonal errors only
	m0dif_header = """
VARNAMES: ROW  zMIN     zMAX     z        MUDIF  MUDIFERR   MUREF  NFIT"""

	lcp = txtobj('lcparam_DS17f.txt') #../cosmosis/cosmosis-standard-library/likelihood/pantheon/
	cov = np.loadtxt('../cosmosis/cosmosis-standard-library/likelihood/pantheon/sys_DS17f.txt',unpack=True,skiprows=1)
	cov = cov.reshape([40,40])

	with open('pantest.m0dif','w') as fout:
		print(m0dif_header,file=fout)
		for i in range(len(lcp.mb)):
			z = lcp.zcmb[i]
			zerr = 150/3e5*5.0/np.log(10)*(1.0+z)/(z*(1.0+z/2.0))
			print(f"ROW:      {i+1}   {lcp.zcmb[i]:.4f} {lcp.zcmb[i]:.4f} {lcp.zcmb[i]:.4f}    0.0000    {np.sqrt(lcp.dmb[i]**2.+cov[i,i]+zerr**2.):.4f}  {lcp.mb[i]+19.36:.4f}    1",file=fout)

if __name__ == "__main__":
	wfit()
