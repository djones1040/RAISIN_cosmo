#!/usr/bin/env python
"plot hubble diagram(s)"

import pylab as plt
plt.ion()
import numpy as np
import os
from txtobj import txtobj
from astropy.cosmology import Planck15
import pdb

def lowz(frfile_opt='output/fit_optical/LOWZ_RAISIN_optical.FITRES.TEXT',
		 frfile='output/fit_nir/LOWZ_RAISIN.FITRES.TEXT',
		 frfile_late='output/fit_nir/LOWZ_RAISIN_latetimes.FITRES.TEXT'):

	fropt = txtobj(frfile_opt,fitresheader=True)
	fr = txtobj(frfile,fitresheader=True)
	frlate = txtobj(frfile_late,fitresheader=True)
	
	plt.clf()
	plt.errorbar(fropt.zCMB,fropt.DLMAG-Planck15.distmod(fropt.zCMB).value,yerr=fropt.DLMAGERR,fmt='o')
	iCut = (fropt.PKMJDERR < 2) & (np.abs(fropt.AV) < 1) #fropt.FITPROB > 0
	#for k in fropt.__dict__.keys():
	#	fropt.__dict__[k] = fropt.__dict__[k][iCut]
	#plt.errorbar(fropt.zCMB,fropt.DLMAG-Planck15.distmod(fropt.zCMB).value,yerr=fropt.DLMAGERR,fmt='o')

	plt.errorbar(fr.zCMB,fr.DLMAG-Planck15.distmod(fr.zCMB).value,yerr=fr.DLMAGERR,fmt='o')

	plt.errorbar(frlate.zCMB,frlate.DLMAG-Planck15.distmod(frlate.zCMB).value,yerr=frlate.DLMAGERR,fmt='o')
	medresid = np.median(frlate.DLMAG-Planck15.distmod(frlate.zCMB).value)
	
	av,stretch,pkmjderr,fitprob = [],[],[],[]
	for j,i in enumerate(frlate.CID):
		if np.abs(frlate.DLMAG[j]-Planck15.distmod(frlate.zCMB[j]).value - medresid) > 0.5:
			try:
				av += [fropt.AV[fropt.CID == i][0]]
				stretch += [fropt.STRETCH[fropt.CID == i][0]]
				pkmjderr += [fropt.PKMJDERR[fropt.CID == i][0]]
				fitprob += [fropt.FITPROB[fropt.CID == i][0]]
			except: continue
				
	plt.axhline(0,lw=2,color='k')
	plt.xlabel('$z_{CMB}$',fontsize=15)
	plt.ylabel('$\mu - \mu_{\Lambda CDM}$',fontsize=15)

	pdb.set_trace()
	
	return

def ps1(frfile_opt='output/fit_optical/PS1_RAISIN_optical.FITRES.TEXT',
		frfile='output/fit_nir/PS1_RAISIN.FITRES.TEXT'):

	fropt = txtobj(frfile_opt,fitresheader=True)
	fr = txtobj(frfile,fitresheader=True)
	
	plt.clf()
	plt.errorbar(fropt.zCMB,fropt.DLMAG-Planck15.distmod(fropt.zCMB).value,yerr=fropt.DLMAGERR,fmt='o')
	iCut = (fropt.PKMJDERR < 2) & (np.abs(fropt.AV) < 1) #fropt.FITPROB > 0
	for k in fropt.__dict__.keys():
		fropt.__dict__[k] = fropt.__dict__[k][iCut]
	plt.errorbar(fropt.zCMB,fropt.DLMAG-Planck15.distmod(fropt.zCMB).value,yerr=fropt.DLMAGERR,fmt='o')

	plt.errorbar(fr.zCMB,fr.DLMAG-Planck15.distmod(fr.zCMB).value,yerr=fr.DLMAGERR,fmt='o')
				
	plt.axhline(0,lw=2,color='k')
	plt.xlabel('$z_{CMB}$',fontsize=15)
	plt.ylabel('$\mu - \mu_{\Lambda CDM}$',fontsize=15)

	pdb.set_trace()
	
	return

def des(frfile_opt='output/fit_optical/DES_RAISIN_optical.FITRES.TEXT',
		frfile='output/fit_nir/DES_RAISIN.FITRES.TEXT'):

	fropt = txtobj(frfile_opt,fitresheader=True)
	fr = txtobj(frfile,fitresheader=True)
	
	plt.clf()
	plt.errorbar(fropt.zCMB,fropt.DLMAG-Planck15.distmod(fropt.zCMB).value,yerr=fropt.DLMAGERR,fmt='o')
	iCut = (fropt.PKMJDERR < 2) & (np.abs(fropt.AV) < 1) #fropt.FITPROB > 0
	for k in fropt.__dict__.keys():
		fropt.__dict__[k] = fropt.__dict__[k][iCut]
	plt.errorbar(fropt.zCMB,fropt.DLMAG-Planck15.distmod(fropt.zCMB).value,yerr=fropt.DLMAGERR,fmt='o')

	plt.errorbar(fr.zCMB,fr.DLMAG-Planck15.distmod(fr.zCMB).value,yerr=fr.DLMAGERR,fmt='o')
				
	plt.axhline(0,lw=2,color='k')
	plt.xlabel('$z_{CMB}$',fontsize=15)
	plt.ylabel('$\mu - \mu_{\Lambda CDM}$',fontsize=15)

	pdb.set_trace()
	
	return


if __name__ == """__main__""":
	#lowz()
	ps1()
	des()
