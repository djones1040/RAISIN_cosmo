#!/usr/bin/env python
# D. Jones - 2/11/2020

import pylab as plt
import numpy as np
from txtobj import txtobj

def mkplots(fr,title=None):
	plt.rcParams['figure.figsize'] = (12,4)
	
	plt.clf()
	ax1,ax2,ax3 = plt.subplot(131),plt.subplot(132),plt.subplot(133)
	dlmagbins = np.linspace(-2,2,25)
	ax1.hist(fr.DLMAG-fr.SIM_DLMAG,bins=dlmagbins)
	ax1.axvline(np.median(fr.DLMAG-fr.SIM_DLMAG),color='r',
				label='median: %.2f'%np.median(fr.DLMAG-fr.SIM_DLMAG))
	stretchbins = np.linspace(-0.2,0.2,25)
	ax2.hist(fr.STRETCH-fr.SIM_STRETCH,bins=stretchbins)
	ax2.axvline(np.median(fr.STRETCH-fr.SIM_STRETCH),color='r',
				label='median: %.2f'%np.median(fr.STRETCH-fr.SIM_STRETCH))
	#import pdb; pdb.set_trace()
	
	avbins = np.linspace(-1,1,25)
	ax3.hist(fr.AV-fr.SIM_AV,bins=avbins)
	ax3.axvline(np.median(fr.AV-fr.SIM_AV),color='r',
				label='median: %.2f'%np.median(fr.AV-fr.SIM_AV))

	ax1.set_xlabel('DLMAG',fontsize=15)
	ax2.set_xlabel('$s_{BV}$',fontsize=15)
	ax3.set_xlabel('$A_V$',fontsize=15)
	if title:
		ax2.set_title(title)
	for ax in [ax1,ax2,ax3]:
		ax.legend()
		
def checksimpar(
		lowzfile='output/sim/LOWZ_RAISIN_SIM.FITRES.TEXT',
		ps1file='output/sim/PS1_RAISIN_SIM.FITRES.TEXT',
		desfile='output/sim/DES_RAISIN_SIM.FITRES.TEXT'):

	frlowz = txtobj(lowzfile,fitresheader=True)
	mkplots(frlowz,title='Low-$z$')
	plt.savefig('figs/checksimfit_lowz.png')
	
	# DES sims
	frdes = txtobj(desfile,fitresheader=True)
	mkplots(frdes,title='DES')
	plt.savefig('figs/checksimfit_des.png')
	
	# PS1 sims
	frps1 = txtobj(ps1file,fitresheader=True)
	mkplots(frps1,title='PS1')
	plt.savefig('figs/checksimfit_ps1.png')
	
if __name__ == "__main__":
	checksimpar()
