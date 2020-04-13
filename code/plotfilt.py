#!/usr/bin/env python

import pylab as plt
import numpy as np
import os
plt.ion()

def main():

	plt.clf()
	ax1 = plt.subplot(131)
	ax2 = plt.subplot(132)
	ax3 = plt.subplot(133)
	
	for filtname,filtfile,low,high,ax in zip(
			['Y','J','H'],
			['$SNDATA_ROOT/filters/CSP/CSP_TAMU/Y.dat',
			 '$SNDATA_ROOT/filters/CSP/CSP_TAMU/J.dat',
			 '$SNDATA_ROOT/filters/CSP/CSP_TAMU/H.dat'],
			[8983,11250,14346],
			[11250,14346,19090],[ax1,ax2,ax3]):
		
		wave,tp = np.loadtxt(os.path.expandvars(filtfile),unpack=True)

		ax.plot(wave,tp)
		ax.axvline(low,color='k')
		ax.axvline(high,color='k')
		ax.set_xlabel('wavelength')
		ax.set_ylabel('%s'%filtname)

	import pdb; pdb.set_trace()

if __name__ == "__main__":
	main()
