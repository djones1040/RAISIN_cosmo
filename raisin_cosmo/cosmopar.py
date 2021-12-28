import getdist.plots as gplot
import pylab as plt
import numpy as np
import pylab as plt
import os

def getw_cosmosis(name='',dirname=None):

	# should run
	# postprocess {inifile} -o plots -p {name} --no-plots
	# if needed
	if dirname is not None:
		cwd = os.getcwd()
		os.chdir(dirname)

	with open(f"plots/{name}_means.txt") as fin:
		for line in fin:
			if line.startswith("cosmological_parameters--omega_m"):
				omegam,omegamerr = float(line.split()[1]),float(line.split()[2])
			elif line.startswith("cosmological_parameters--w"):
				w,werr = float(line.split()[1]),float(line.split()[2])

	if dirname is not None:
		os.chdir(cwd)

	return w,werr,omegam,omegamerr

#def getw_cosmosis(name=''):

	# should run
	# postprocess {inifile} -o plots -p {name} --no-plots
	# if needed
#	 with open(f"plots/{name}_means.txt") as fin:
#		 for line in fin:
#			 if line.startswith("cosmological_parameters--omega_m"):
#				 omegam,omegamerr = float(line.split()[1]),float(line.split()[2])
#			 elif line.startswith("cosmological_parameters--w"):
#				 w,werr = float(line.split()[1]),float(line.split()[2])
#	 return omegam,omegamerr,w,werr

def getw(name=''):

	g = gplot.getSinglePlotter(chain_dir='/scratch/midway2/rkessler/djones/cosmomc/chains_2015/chains/')
	samples = g.sampleAnalyser.samplesForRoot(name)

	p = samples.getParams()

	plt.close()
	print(f"w: {samples.mean(p.w)} +/- {samples.std(p.w)}")
	print(f"Om: {samples.mean(p.omegam)} +/- {samples.std(p.omegam)}")
	return(samples.mean(p.w),samples.std(p.w),
		   samples.mean(p.omegam),samples.std(p.omegam))

def geth0(name=''):

	g = gplot.getSinglePlotter(chain_dir='/scratch/midway2/rkessler/djones/cosmomc/chains_2015/chains/')
	samples = g.sampleAnalyser.samplesForRoot(name)

	p = samples.getParams()

	plt.close()
#	 import pdb; pdb.set_trace()
	print(f"H0: {samples.mean(p.__dict__['H0'])} +/- {samples.std(p.__dict__['H0'])}")
	print(f"w: {samples.mean(p.w)} +/- {samples.std(p.w)}")
	print(f"Om: {samples.mean(p.omegam)} +/- {samples.std(p.omegam)}")
	return(samples.mean(p.w),samples.std(p.w),
		   samples.mean(p.omegam),samples.std(p.omegam))


def getom(name=''):

	#g = gplot.getSinglePlotter(chain_dir='/scratch/midway2/rkessler/djones/cosmomc/chains_2015/')
	g = gplot.getSinglePlotter(chain_dir='/n/holystore01/LABS/berger_lab/Lab/djones01/RAISIN/cosmo/cosmomc_chains')
	samples = g.sampleAnalyser.samplesForRoot(name)

	p = samples.getParams()

	plt.close()
	print(f"Om: {samples.mean(p.omegam)} +/- {samples.std(p.omegam)}")
	return(samples.mean(p.omegam),samples.std(p.omegam))

def getom_odyssey(name=''):

	g = gplot.getSinglePlotter(chain_dir='/n/holystore01/LABS/berger_lab/Lab/djones01/RAISIN/cosmo/cosmomc/chains_odyssey')
	samples = g.sampleAnalyser.samplesForRoot(name)

	p = samples.getParams()

	plt.close()
	print(f"Om: {samples.mean(p.omegam)} +/- {samples.std(p.omegam)}")
	return(samples.mean(p.omegam),samples.std(p.omegam))


def cosmosys(postprocess=False,version='nir'):

	syslist = ['stat','all','photcal','massdivide','biascor','pecvel','kcor','mwebv','lowzcal','hstcal','tmpl','lcfitter']
	titles = ['Stat.','All Sys.','Phot. Cal.','Mass Step',
			  'Bias Corr.','Pec. Vel.','$k$-corr.',
			  'MW E(B-V)','Low-$z$ Cal','$HST$ Cal','Template Flux','NIR SN Model']
	fileprefix='raisin'


	if postprocess:
		os.system("postprocess RAISIN_stat.ini -o plots -p raisin_stat --no-plots")
		os.system("postprocess RAISIN_all.ini -o plots -p raisin_all --no-plots")

		for s in syslist:
			os.system(f"postprocess RAISIN_{s}.ini -o plots -p raisin_{s} --no-plots")
		
	wstat,werrstat,omstat,omerrstat = getw_cosmosis(name='raisin_stat')
	wall,werrall,omall,omerrall = getw_cosmosis(name='raisin_all')

	tblhdr = """\\begin{deluxetable}{lccc}
\\tabletypesize{\\scriptsize}
\\tablewidth{0pt}
\\tablecaption{Summary of Systematic Uncertainties on $w$}
\\tablehead{Error&$\\Delta w$&$\\sigma_w$&$\\sigma_w/\\sigma_{\\textrm{stat}}$}"""
	print(tblhdr)

	for s,t,i in zip(syslist,titles,range(len(titles))):
		f = '%s_%s'%(fileprefix,s)
		w,werr,om,omerr = getw_cosmosis(name=f)
		if werr < werrstat: 
			werrsys = 0.0
		else:
			werrsys = np.sqrt(werr**2. - werrstat**2.)
		if f == 'foundps1_all':
			print('%s&%.3f&%.3f&%.3f\\\\'%(
				t,wall-wstat,werrsys,werrsys/werrstat))
		else:
			print('%s&%.3f&%.3f&%.3f\\\\'%(
				t,w-wstat,werrsys,werrsys/werrstat))
#		 import pdb; pdb.set_trace()
		if i == 0:
			print('\\tableline\\\\')

		tblfooter = """\\startdata
\\enddata
\\label{table:syserr}
\\end{deluxetable}"""
	print(tblfooter)

def syspiechart(ax=None,
				sysval=[0.025,0.075,0.039,0.012,0.011],
				title=None,
				syslist=['Phot. Cal.','Bias Corr.', #'$k$-corr.',
						 #'Pec. Vel.',
						 'Mass\nStep',#'NIR\nModel',
						 'Template\nFlux',
						 'Other'],
				explode=[0,0,0,0,0,0],radius=1.4,fontsize=13,makebold=False,startangle=65):
	import matplotlib.patheffects as path_effects

	#			 sysval=[0.027,0.044,0.005,0.020,0.043,0.004,0.004,0.002],
	#			 title=None,
	#			 syslist=['Phot. Cal.','Bias Corr.','$k$-corr.',
	#					  'Pec. Vel.','Mass\nStep','NIR SN\nModel','MW E(B-V)','Template Flux'],

	
	systot = np.sum(sysval)
	sysval /= np.sum(sysval)
	colors = ['#fbb4ae',
			  '#b3cde3',
			  '#ccebc5',
			  '#decbe4',
			  '#fed9a6',
			  '#ffffcc',
			  '#e5d8bd',
			  '#fddaec']

	def absolute_value(val):
		#a	  = np.round(val/np.sum(sysval), 0)
		if val*systot/100 > 0.004:
			return '%.3f\n(%i%%)'%(val*systot/100,val)
		else:
			return ''
	plt.rcParams['font.size'] = 17
	if not ax:
		plt.clf()
		ax = plt.axes()
	patches, texts, autotexts = ax.pie(sysval, labels=syslist, colors=colors,
									   autopct=absolute_value, shadow=False, startangle=startangle,
									   labeldistance=1.08,#explode=explode,
									   #wedgeprops = {'linewidth': 2, 'edgecolor':'k'},
									   pctdistance=0.55,radius=radius)
	ax.set_title(title,y=1.025)
	for patch in patches:
		patch.set_path_effects([path_effects.SimpleLineShadow(),
								path_effects.Normal()])
		
	for label,text in zip(autotexts,texts):
		if makebold:
			if 'CC' in text._text or 'SALT2' in text._text or 'Bias' in text._text:
				label._fontproperties._weight = 'bold'
				text._fontproperties._weight = 'bold'
		label._fontproperties._size = fontsize
		#import pdb; pdb.set_trace()

	centre_circle = plt.Circle((0,0),1.00,fc='white')
	fig = plt.gcf()
	fig.gca().add_artist(centre_circle)
	plt.savefig('raisin_syspie.png',dpi=200)

def getcorner(name=''):

	import corner

	g = gplot.getSinglePlotter(chain_dir='/scratch/midway2/rkessler/djones/cosmomc/chains_2015/chains/')
	samples = g.sampleAnalyser.samplesForRoot(name)

	p = samples.getParams()

	mat = np.zeros([len(p.w),3])
	mat[:,0] = p.__dict__['H0']
	mat[:,1] = p.w
	mat[:,2] = p.omegam
	corner.corner(mat,labels=['H$_0$','$w$','$\Omega_m$'],
				  quantiles=[0.16,0.5,0.84],show_titles=True,
				  title_kwargs={"fontsize":20},label_kwargs={"fontsize":20},bins=40)
	plt.savefig('cosmo_corner.png',dpi=200)
	#plt.close()

def getcorner_cosmosis(name=''):

    import corner

	#g = gplot.getSinglePlotter(chain_dir='/scratch/midway2/rkessler/djones/cosmomc/chains_2015/chains/')
	#samples = g.sampleAnalyser.samplesForRoot(name)

	#p = samples.getParams()
	omegam,w,h0,weights = np.loadtxt('raisin_cosmosis_chains.txt',unpack=True)
	weights /= weights.max()
	h0 *= 100
	
	iplot = (h0 > 62) & (h0 < 75) & (w > -1.33) & (w < -0.83) & (omegam > 0.23) & (omegam < 0.4)
	omegam,w,h0,weights = omegam[iplot],w[iplot],h0[iplot],weights[iplot]
	#import pdb; pdb.set_trace()
	
	mat = np.zeros([len(w),3])
	mat[:,0] = h0
	mat[:,1] = w
	mat[:,2] = omegam
	corner.corner(mat,labels=['H$_0$','$w$','$\Omega_m$'],
				  quantiles=[0.16,0.5,0.84],show_titles=True,weights=weights,
				  title_kwargs={"fontsize":20},label_kwargs={"fontsize":20},bins=30,plot_datapoints=False,
				  plot_density=True,smooth=1.0,no_fill_contours=True,fill_contours=False)

	#fig = corner.corner(mat, plot_datapoints=True, labels=['H$_0$','$w$','$\Omega_m$'],
	#					 fill_contours=False, weights=weights, levels=[0.68, 0.95], bins=25,
	#					 smooth=1.0, no_fill_contours=True, plot_density=False,
	#					 color='k',show_titles='True')
	plt.savefig('cosmo_corner.png',dpi=200)
	#plt.close()

def cosmotable_optnir():
	# comparing cosmological results from optical/NIR variants
	for var,label in zip(
			['optical','nir','opticalnir'],
			['Optical only','NIR only','Optical$+$NIR']):
		
		# get the cosmo params
		wstat,werrstat,omstat,omerrstat = getw_cosmosis(name='raisin_stat',dirname=f'cosmosis_{var}')
		wall,werrall,omall,omerrall = getw_cosmosis(name='raisin_all',dirname=f'cosmosis_{var}')
		werrsys = np.sqrt(werrall**2. - werrstat**2.)
		print(f"{label}&${wstat:.3f}\\pm{werrstat:.3f}$&${wall:.3f}\\pm{werrall:.3f}$&${werrsys/werrstat:.3f}$\\\\")

if __name__ == "__main__":

	#getw('planck_pan_wcdm_approx')
	#getw('raisin_all_planck18')
	#getw('planck18')
	#getw('raisin_wcdm_snalone')

    #getw('raisin_stat')
    #geth0('raisin_all')
    getcorner_cosmosis('raisin_all')
    #getom('RAISIN_all_ocdm')
    #getom('RAISIN_all_lcdm')
    getom('RAISIN_combined_all_lcdm')
	#getom('planck_lcdm_approx')
	#getw('sn_cmb_omw_0')
	#getom_odyssey('RAISIN_combined_all_lcdm')

	#cosmosys(postprocess=False)
	#om,omerr,w,werr = getw_cosmosis('raisin_stat')
	#print(w,werr)
	#w,werr,om,omerr = getw_cosmosis('raisin_all')
	#print(w,werr)
	
	#syspiechart()
	#cosmotable_optnir()

