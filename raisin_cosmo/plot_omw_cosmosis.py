import getdist.plots as gplot
import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import numpy as np
import plotsetup
plotsetup.halfpaperfig()

g = gplot.getPlotter(chain_dir = '/Users/david/Dropbox/research/raisin/cosmo/cosmomc_chains')
#/scratch/midway/rkessler/dscolnic/cosmomc')
g.settings.fig_width_inch = 5.0
g.settings.lw_contour = 2.0
g.settings.solid_contour_palefactor = 0.2
g.settings.legend_frame = False
g.settings.axes_fontsize=12
g.settings.legend_fontsize=12
g.settings.lab_fontsize=14  
#roots = ['chains/noSN_omw', 'chains2/DS17_ALL_omw_alone','chains2/DS17_ALL_nosys_omw_alone','chains/noSN_bao_omw','chains2/DS17_ALL_omw','chains2/DS17_ALL_nosys_omw']
#names=['CMB','SN','SN(stat)','CMB+BAO','SN+CMB','SN(stat)+CMB']  
#colors=['teal','k','gray','blue','purple','gray']


roots = ['raisin_wcdm_snalone','DS17_ALL_omw_alone']
names=['CMB','SN']
colors=['darkcyan','#C1292E']#,'C1','#020100','#235789'] 

g.plot_2d(roots,'omegam','w',
		  filled=[True,True],colors=colors)

plt.xlabel(r'$\Omega_m$',labelpad=10)
plt.ylabel(r'$w$',labelpad=5)
plt.xlim(0.15,0.475)
plt.ylim(-1.5,-.6)
plt.title(r'$w \rm{CDM}$'+' Constraints For Combined Samples')

plt.text(0.3,-1.45,"$Planck$+18\nalone", color='b',rotation=0,ha='center',va='center',
		 bbox={'facecolor':'1.0','edgecolor':'1.0','boxstyle':'round','alpha':0.5})
plt.text(0.257,-0.9,"Pantheon\nSNe", color='k',rotation=0,alpha=0.5,ha='center',va='center',
		 bbox={'facecolor':'1.0','edgecolor':'1.0','boxstyle':'round','alpha':0.5})
#plt.text(0.29,-1.23,"Pantheon SN (Stat)", color='gray',rotation=305,alpha=0.5)
#plt.text(0.28,-.75,"CMB+BAO", color='yellow')
props = dict(boxstyle='round', facecolor='yellow', alpha=0.5)
#plt.text(0.26,-.8,"CMB+BAO", color='blue', bbox=props)

props = dict(boxstyle='round', facecolor='black', alpha=0.5)
plt.text(0.417,-1.06,"Pantheon SNe \n+ CMB",color='white', bbox=props,ha='center',va='center')

props = dict(boxstyle='round', facecolor='C1', alpha=0.8)
plt.text(0.212,-1.2,"RAISIN +\n CMB",color='k', bbox=props,ha='center',va='center')

props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
plt.text(0.205,-1.075,"RAISIN SNe\nalone",color='k', bbox=props,ha='center',va='center')

ax = plt.gca()
str2=colors

gcf = plt.figure(figsize=(20,10))
str2=colors

# ok, now let's try the cosmosis stuff.....
#from cosmosis.postprocessing.postprocess import postprocessor_for_sampler
from cosmosis.postprocessing.inputs import read_input
import cosmosis.postprocessing.postprocess
import cosmosis.postprocessing.plots

#sampler, ini = read_input(ini_filename, args.text, args.weights)
#processor_class = postprocessor_for_sampler(sampler.split()[-1])
#processor = processor_class(ini, 'whatever', i) #, **vars(args))
#processor.run()
#processor.finalize()

sampler, ini = read_input('cosmosis/planck_samples/chain_p-TTTEEE-lowE_wcdm.ini')
ip = cosmosis.postprocessing.postprocess.ImportanceProcessor(ini=ini,label='raisin',index=0)
mp = cosmosis.postprocessing.plots.WeightedMetropolisPlots2D(data_source=ip)
mp.plot_set = 0
mp.make_2d_plot('cosmological_parameters--omega_m','cosmological_parameters--w',ax=ax)

sampler, ini = read_input('cosmosis/RAISIN_all.ini')
ip = cosmosis.postprocessing.postprocess.ImportanceProcessor(ini=ini,label='raisin',index=0)
mp = cosmosis.postprocessing.plots.WeightedMetropolisPlots2D(data_source=ip)
mp.plot_set = 1
mp.make_2d_plot('cosmological_parameters--omega_m','cosmological_parameters--w',ax=ax,alpha=0.75)
#mp.make_2d_plot('cosmological_parameters--h0','cosmological_parameters--w',ax=ax,alpha=0.75)

#sampler, ini = read_input('cosmosis/pan_sne.ini')
#ip = cosmosis.postprocessing.postprocess.ImportanceProcessor(ini=ini,label='raisin',index=0)
#mp = cosmosis.postprocessing.plots.WeightedMetropolisPlots2D(data_source=ip)
#mp.plot_set=2
#mp.make_2d_plot('cosmological_parameters--omega_m','cosmological_parameters--w',ax=ax,alpha=0.5)

ax.axhline(-1.0,color='k',lw=2)


##########
g.export('plot_omw.png')

