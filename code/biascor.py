#!/usr/bin/env python
# D. Jones - 2/27/20

import numpy as np
from scipy.stats import binned_statistic
from txtobj import txtobj
import pylab as plt
plt.ion()
# rsync -avz output/sim/PS1_RAISIN_SIM_NIR david@turanga.ucsc.edu:/Users/David/Dropbox/research/RAISIN/cosmo/output/sim/

def errfnc(x):
    return(np.std(x)/np.sqrt(len(x)))

def main(frsim_ps1='output/sim/PS1_RAISIN_SIM.FITRES.TEXT',#'output/sim/PS1_RAISIN_SIM_NIR.FITRES.TEXT',
		 #output/sim/PS1_RAISIN_SIM_NIR/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
		 frsim_des='output/sim/DES_RAISIN_SIM.FITRES.TEXT'): #'output/sim/DES_RAISIN_SIM_NIR.FITRES.TEXT'):

	plt.clf()
	ax = plt.axes() #subplot(121)
	
	frp = txtobj(frsim_ps1,fitresheader=True)
	frd = txtobj(frsim_des,fitresheader=True)

	zbins = np.linspace(0.01,0.7,20)
	#iSNR = frp.SNRMAX1 > 50
	#ps1_bin = binned_statistic(
	#	frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
	#	statistic='median',bins=zbins).statistic
	#ps1_errbin = binned_statistic(
	#	frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
	#	statistic=errfnc,bins=zbins).statistic
	#ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
	#			fmt='o-',color='b',capsize=0,lw=2)

	ps1_bin = binned_statistic(
		frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
		statistic='median',bins=zbins).statistic
	ps1_errbin = binned_statistic(
		frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
		statistic=errfnc,bins=zbins).statistic

	
	des_bin = binned_statistic(
		frd.zCMB,frd.DLMAG-frd.SIM_DLMAG,
		statistic='median',bins=zbins).statistic

	#ax2 = plt.subplot(122)
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
				fmt='o-',color='r',capsize=0,lw=2)
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,des_bin,
				fmt='o-',color='b',capsize=0,lw=2)
	ax.set_xlabel('$s_{CMB}$',fontsize=15)
	ax.set_ylabel('$\mu - \mu_{sim}$',fontsize=15)
	import pdb; pdb.set_trace()
	
	return
	
if __name__ == "__main__":
	main()
