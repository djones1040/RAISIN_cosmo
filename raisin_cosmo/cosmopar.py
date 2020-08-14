import getdist.plots as gplot
import pylab as plt
import numpy as np
import pylab as plt

def getw(name=''):

	g = gplot.getSinglePlotter(chain_dir='/scratch/midway2/rkessler/djones/cosmomc/chains/')
	samples = g.sampleAnalyser.samplesForRoot(name)

	p = samples.getParams()

	plt.close()
	print(f"w: {samples.mean(p.w)} +/- {samples.std(p.w)}")
	print(f"Om: {samples.mean(p.omegam)} +/- {samples.std(p.omegam)}")
	return(samples.mean(p.w),samples.std(p.w),
		   samples.mean(p.omegam),samples.std(p.omegam))

def getom(name=''):

	g = gplot.getSinglePlotter(chain_dir='/scratch/midway2/rkessler/djones/cosmomc/chains/')
	samples = g.sampleAnalyser.samplesForRoot(name)

	p = samples.getParams()

	plt.close()
	print(f"Om: {samples.mean(p.omegam)} +/- {samples.std(p.omegam)}")
	return(samples.mean(p.omegam),samples.std(p.omegam))


if __name__ == "__main__":
	#getw('RAISIN_stat')
	#getom('RAISIN_stat_lcdm')
	getw('sn_cmb_omw_0')

