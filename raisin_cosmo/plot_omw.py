import getdist.plots as gplot
import matplotlib.pyplot as plt
import numpy as np
import plotsetup
plotsetup.halfpaperfig()

g = gplot.getPlotter(chain_dir = '/scratch/midway2/rkessler/dscolnic/cosmomc/')
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


roots = ['chains2/noSN_omw','chains2/DS17_ALL_omw_alone','chains2/DS17_ALL_omw','chains_raisin/RAISIN_stat']
names=['CMB','SN','SN+CMB','RAISIN']
colors=['teal','k','purple']
colors=['#235789','#C1292E','#020100','C1'] 

#roots=['chains2/DS17_ALL_omw','chains/noSN_bao_omw','chains2/DS17_ALL_nosys_omw_alone','chains2/DS17_ALL_omw_alone','chains/noSN_omw']
#colors=['#020100','#F1D302','#FDFFFC','#C1292E','#235789']

#g.add_y_bands(-1,0.01,alpha1=0.15)
#ax = plt.axes()
g.plot_2d(roots,'omegam','w',
		  filled=[True,True,True,True],colors=colors)
#g2.plot_2d('chains_raisin/RAISIN_stat','omegam','w',
#		   filled=[True],colors=['r'],subplots=ax)

#g.add_y_marker(-1)         
#g.plot_2d(roots, 'omegam', 'w', filled=True)
#g.legend(['CMB', 'CMB+BAO','SN+CMB','SN'], legend_loc='upper left',frameon=False,fontsize=8,ncol=2);
#line, = plt.plot(np.arange(-5,5,0.1), np.arange(-5,5,0.1)*0-1, lw=1,alpha=0.1)

plt.xlabel(r'$\Omega_m$',labelpad=10)
plt.ylabel(r'$w$',labelpad=5)
plt.xlim(0.15,0.5)
plt.ylim(-1.5,-.5)
plt.title(r'$w \rm{CDM}$'+' Constraints For Combined Samples')

plt.text(0.42,-.68,"CMB", color='blue',rotation=40,ha='center',va='center',
		 bbox={'facecolor':'1.0','edgecolor':'1.0','boxstyle':'round','alpha':0.5})
plt.text(0.37,-1.25,"Pantheon SN", color='red',rotation=306,alpha=0.5,ha='center',va='center',
		 bbox={'facecolor':'1.0','edgecolor':'1.0','boxstyle':'round','alpha':0.5})
#plt.text(0.29,-1.23,"Pantheon SN (Stat)", color='gray',rotation=305,alpha=0.5)
#plt.text(0.28,-.75,"CMB+BAO", color='yellow')
props = dict(boxstyle='round', facecolor='yellow', alpha=0.5)
#plt.text(0.26,-.8,"CMB+BAO", color='blue', bbox=props)

props = dict(boxstyle='round', facecolor='black', alpha=0.8)
plt.text(0.35,-1.0,"Pantheon SN + CMB",color='white', bbox=props)

props = dict(boxstyle='round', facecolor='C1', alpha=0.8)
plt.text(0.22,-1.2,"RAISIN (stat) +\n CMB",color='k', bbox=props,ha='center',va='center')


#plt.xticks([])
#plt.yticks([])
str2=colors
#colors=['black','green','blue','red','purple']
#proxy = [plt.Rectangle((0,0),.05,.05,fc = str2[pc]) for pc in range(0,3)]
#names=['CMB','SN','SN(stat)','CMB+BAO','SN+CMB','SN(stat)+CMB']
#plt.legend(proxy,names,loc='upper left',prop={'size':12},frameon=False,ncol=3)

plt.figure(figsize=(20,10))
str2=colors
#colors=['black','green','blue','red','purple']
#'CMB', 'CMB+BAO','SN+CMB','SN'


g.export('plot_omw.png')
#plt.savefig('plot_omw.png')
