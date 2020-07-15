#!/usr/bin/env python
# D. Jones - 7/14/20

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
plt.ion()
import numpy as np
from txtobj import txtobj
from scipy.optimize import minimize
import scipy.stats
import os
import cosmo
plt.rcParams['figure.figsize'] = (11,6)
#plt.subplots_adjust(hspace=0)

_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_CIDS.LIST',unpack=True,dtype=str),
							#np.loadtxt('output/goodcids/CfA_CIDS.LIST',unpack=True,dtype=str),
							np.loadtxt('output/goodcids/PS1_CIDS.LIST',unpack=True,dtype=str),
							np.loadtxt('output/goodcids/DES_CIDS.LIST',unpack=True,dtype=str)))

nirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN.FITRES.TEXT',
					 '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN.FITRES.TEXT',
					 '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN.FITRES.TEXT']
nirdatafitreslist = [os.path.expandvars(filepath) for filepath in nirdatafitreslist]
nirdatafitreslist = ['output/fitres_cosmo/CSP.FITRES',
					 'output/fitres_cosmo/PS1.FITRES',
					 'output/fitres_cosmo/DES.FITRES']

opticalnirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',
							'$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT',
							'$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT']
opticalnirdatafitreslist = [os.path.expandvars(filepath) for filepath in opticalnirdatafitreslist]

def format_axes(fig):
	for i, ax in enumerate(fig.axes):
		ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
		ax.tick_params(labelbottom=False, labelleft=False)
		
def lnlikefunc(x,p_iae=None,mu_i=None,sigma_i=None,sigma=None):

	if sigma or sigma == 0.0:
		# fix the dispersion
		x[2] = sigma; x[3] = sigma

	p_iae[np.where(p_iae == 0)] == 1e-4
	return -np.sum(np.logaddexp(-(mu_i-x[0])**2./(2.0*(sigma_i**2.+x[2]**2.)) +\
				np.log((1-0.01*p_iae)/(np.sqrt(2*np.pi)*np.sqrt(x[2]**2.+sigma_i**2.))),
				-(mu_i-x[1])**2./(2.0*(sigma_i**2.+x[3]**2.)) +\
				np.log((0.01*p_iae)/(np.sqrt(2*np.pi)*np.sqrt(x[3]**2.+sigma_i**2.)))))

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

	#iGoodSt = (fropt.STRETCH > 0.8) & (fropt.STRETCH < 1.3)

	for k in fr.__dict__.keys():
		fr.__dict__[k] = fr.__dict__[k][iGoodAV & iGoodSt]

	#restrict_to_good_list = False
	if restrict_to_good_list:
		# then case-by-case removal of bad stuff
		# because some RAISIN things have bad photometry
		iGood = np.array([],dtype=int)
		for j,i in enumerate(fr.CID):
			if i in _goodcids: iGood = np.append(iGood,j)
		for k in fr.__dict__.keys():
			fr.__dict__[k] = fr.__dict__[k][iGood]

	return fr


def main(boundary=10):

	fig = plt.figure()#constrained_layout=True)
	plt.subplots_adjust(top=0.95)
	gs = GridSpec(3, 3, figure=fig)
	gs.update(wspace=0.0, hspace=0.4)

	#plt.tick_params(left)
	
	axmain = fig.add_subplot(gs[1:, :])
	ax1 = fig.add_subplot(gs[0, 0])
	ax2 = fig.add_subplot(gs[0, 1])
	ax3 = fig.add_subplot(gs[0, 2])


	mp_full,mass_full,masserr_full,resid_full,residerr_full = \
		np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
	
	for frfile,froptfile,ax,title in zip(
			nirdatafitreslist,opticalnirdatafitreslist,[ax1,ax2,ax3],['CSP','PS1','DES']):

		fr = txtobj(frfile,fitresheader=True)
		fropt = txtobj(froptfile,fitresheader=True)
		fr = apply_all_cuts(fr,fropt,restrict_to_good_list=True)
		fr.resid = fr.DLMAG - cosmo.mu(fr.zCMB)
		iGood = np.where(fr.HOST_LOGMASS_ERR < 5)[0]
		for k in fr.__dict__.keys():
			fr.__dict__[k] = fr.__dict__[k][iGood]
		
		fr.p_hm = np.zeros(len(fr.CID))
		for i in range(len(fr.CID)):
			fr.p_hm[i] = scipy.stats.norm.cdf(
				boundary,float(fr.HOST_LOGMASS[i]),
				float(fr.HOST_LOGMASS_ERR[i]))*100.
		
		md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
					  args=(fr.p_hm,fr.resid,fr.DLMAGERR,None))

		resid_iaa,resid_iae = md.x[0],md.x[1]
		scat_iaa,scat_iae = md.x[2],md.x[3]
		residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
		covar = np.sqrt(np.abs(md.hess_inv[1,0]))
		step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

		ax.plot(np.arange(boundary-10,boundary,0.001),
				[resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
				lw=2,color='0.2')
		ax.plot(np.arange(boundary,boundary+10,0.001),
				[resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
				boundary,10,lw=2,color='0.2')
		
		ax.fill_between(np.arange(boundary-10,boundary,0.001),
						[resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
						[resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
		ax.fill_between(np.arange(boundary,boundary+10,0.001),
						[resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
						[resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
						color='red',alpha=0.4)

		# light shading for the dispersion
		ax.fill_between(np.arange(boundary-10,boundary,0.001),
						[resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
						[resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
						color='lightblue',zorder=1,alpha=0.4)
		ax.fill_between(np.arange(boundary,boundary+10,0.001),
						[resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
						[resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
						color='red',zorder=1,alpha=0.6)

		ax.axvline(boundary,ls='--',lw=2,color='0.2')
		ax.errorbar(fr.HOST_LOGMASS,fr.resid,xerr=fr.HOST_LOGMASS_ERR,
					yerr=fr.DLMAGERR,color='0.6',fmt='',ls='None')
		sc = ax.scatter(fr.HOST_LOGMASS,fr.resid,c=100-fr.p_hm,
						s=30,zorder=9,cmap='RdBu_r')
		ax.text(0.05,0.95,f"$\Delta_M$ = {step:.2f}$\pm${steperr:.2f}",
				va='top',ha='left',transform=ax.transAxes,bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
				zorder=100)

		mp_full = np.append(mp_full,fr.p_hm)
		mass_full = np.append(mass_full,fr.HOST_LOGMASS)
		masserr_full = np.append(masserr_full,fr.HOST_LOGMASS_ERR)
		resid_full = np.append(resid_full,fr.resid)
		residerr_full = np.append(residerr_full,fr.DLMAGERR)
		
		ax.set_ylim([-0.5,0.5])
		ax.set_xlim([7,13])
		ax.set_title(title)
		ax.set_xlabel('log(M/M$_{\odot}$)')
		ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
		ax.xaxis.set_ticks([8,9,10,11,12])
	ax1.set_ylabel('Hubble Resid')
	ax2.yaxis.set_ticklabels([])
	ax3.tick_params(top="on",bottom="on",left="off",right="on",direction="inout",length=8, width=1.5)
	ax3.yaxis.tick_right()
	


	md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
				  args=(mp_full,resid_full,residerr_full,None))

	resid_iaa,resid_iae = md.x[0],md.x[1]
	scat_iaa,scat_iae = md.x[2],md.x[3]
	residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
	covar = np.sqrt(np.abs(md.hess_inv[1,0]))
	step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

	axmain.plot(np.arange(boundary-10,boundary,0.001),
				[resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
				lw=2,color='0.2')
	axmain.plot(np.arange(boundary,boundary+10,0.001),
				[resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
				boundary,10,lw=2,color='0.2')
		
	axmain.fill_between(np.arange(boundary-10,boundary,0.001),
						[resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
						[resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
	axmain.fill_between(np.arange(boundary,boundary+10,0.001),
						[resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
						[resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
						color='red',alpha=0.4)

	# light shading for the dispersion
	axmain.fill_between(np.arange(boundary-10,boundary,0.001),
						[resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
						[resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
						color='lightblue',zorder=1,alpha=0.4)
	axmain.fill_between(np.arange(boundary,boundary+10,0.001),
						[resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
						[resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
						color='red',zorder=1,alpha=0.6)

	axmain.axvline(boundary,ls='--',lw=2,color='0.2')
	axmain.errorbar(mass_full,resid_full,xerr=masserr_full,
					yerr=residerr_full,color='0.6',fmt='',ls='None')
	sc = axmain.scatter(mass_full,resid_full,c=100-mp_full,
					s=30,zorder=9,cmap='RdBu_r')
	axmain.text(0.02,0.95,f"Total $\Delta_M$ = {step:.2f}$\pm${steperr:.2f}",
				va='top',ha='left',transform=axmain.transAxes,
				bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
				zorder=100,fontsize=15)

	axmain.set_ylim([-0.5,0.5])
	axmain.set_xlim([7,13])
	axmain.set_ylabel('Hubble Resid',fontsize=15)
	axmain.set_xlabel('log(M/M$_{\odot}$)',fontsize=15)
	axmain.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
	axmain.xaxis.set_ticks([8,9,10,11,12])
	plt.savefig('figs/raisin_massstep.png',dpi=200)
	
		
	import pdb; pdb.set_trace()
	
if __name__ == "__main__":
	main()
