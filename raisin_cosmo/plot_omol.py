import getdist.plots as gplot
import matplotlib.pyplot as plt
import numpy as np
import plotsetup

bigmat=np.zeros((65,100))
sys_values=open('mlcs_cont.txt','r').readlines()
for x in range(0,len(sys_values)):
               y=sys_values[x].split()
               print(y)
               #stop
               for z in range(0,len(y)):
                              bigmat[x,z]=float(y[z])

#bigmat=np.zeros((65,100))
om=np.arange(0,2.6,0.04)
ol=np.arange(-1,3.0,0.04)
                                                                                                         
#mlcs_cont_ol.txt
# mlcs_cont_om.txt
#  mlcs_cont.txt

  
plotsetup.halfpaperfig()
print(len(om))
print(len(ol))
#stop
print(om)
print(ol)
z= np.matrix.transpose(bigmat)

#CS = plt.contour(om,ol,bigmat.transpose)
#levels=[2.31, 6.17, 11.8]
g = gplot.getPlotter(chain_dir = '/scratch/midway2/rkessler/dscolnic/cosmomc/')
#'chains2x/scratch/midway/rkessler/dscolnic/cosmomc')
g.settings.fig_width_inch = 5.0
g.settings.lw_contour = 2.0
g.settings.solid_contour_palefactor = 0.2
g.settings.legend_frame = False
g.settings.axes_fontsize=12
g.settings.legend_fontsize=12
g.settings.lab_fontsize=14

#plt.savefig('plot_omol.png')

#stop
#roots = ['chains/noSN_omw', 'chains2/DS17_ALL_omw_alone','chains2/DS17_ALL_nosys_omw_alone','chains/noSN_bao_omw','chains2/DS17_ALL_omw','chains2/DS17_ALL_nosys_omw']
#names=['CMB','SN','SN(stat)','CMB+BAO','SN+CMB','SN(stat)+CMB']  
#colors=['teal','k','gray','blue','purple','gray']

roots = ['chains2/noSN_omol', 'chains2/DS17_ALL_omol_alone','chains2/DS17_ALL_omol']
#roots = ['chains2/noSN_omol','chains2/noSN_omol','chains2/noSN_omol','chains2/noSN_omol']
names=['CMB','SN','SN+CMB']
#names=['SN','SN(stat)','SN+CMB']

roots=['chains_raisin/RAISIN_stat_ocdm','chains2/DS17_ALL_omol_alone']
colors=['#C1292E','RAISIN (stat)']
#colors=['teal','k','blue','purple']
#'#778899'
#colors=['#235789','#C1292E','#778899','#020100','#F1D302','#FDFFFC']
colors=['C1','red']
x=np.arange(0,2,.1)
g.plot_2d(roots,'omegam','omegal*',
          filled=[True,True,False,True],colors=colors)
line, = plt.plot(x, 1-x, lw=2,color='red',linestyle='--',alpha=0.3)
line, = plt.plot(x, 0.5*x, lw=2,color='blue',linestyle='--',alpha=0.3)
#g.plot_2d(roots, 'omegam', 'w', filled=True)
#g.legend(['CMB', 'CMB+BAO','SN+CMB','SN'], legend_loc='upper left',frameon=False,fontsize=8,ncol=2);
#line, = plt.plot(np.arange(-5,5,0.1), np.arange(-5,5,0.1)*0-1, lw=1,alpha=0.1)

plt.xlabel(r'$\Omega_m$',labelpad=10)
plt.ylabel(r'$\Omega_{\Lambda}$',labelpad=5)
plt.xlim(0.0,1.6)
plt.ylim(0.0,2.4)
plt.title(r'$o \rm{CDM}$'+' Constraints For SN-only Sample')

#plt.xticks([])
#plt.yticks([])
str2=colors
#colors=['black','green','blue','red','purple']
#proxy = [plt.Rectangle((0,0),.05,.05,fc = str2[pc]) for pc in range(0,4)]
proxy = [plt.Rectangle((0,0),.05,.05,fc = str2[pc]) for pc in range(0,2)] 
#names=['CMB','SN','SN(stat)','CMB+BAO','SN+CMB','SN(stat)+CMB']


#plt.text(0.13,0.29,"Accelerating Universe", color='blue',rotation=28,alpha=0.5)
#plt.text(0.17,0.23,"Decelerating Universe", color='blue',rotation=28,alpha=0.5)

#plt.text(0.7,0.21,"Flat Universe", color='red',rotation=315,alpha=0.5)

plt.text(1.2,0.67,"Accelerating Universe", color='blue',rotation=18,alpha=0.5,ha='center',va='center')
plt.text(1.22,0.53,"Decelerating Universe", color='blue',rotation=18,alpha=0.5,ha='center',va='center')

plt.text(0.65,0.24,"Flat Universe", color='red',rotation=325,alpha=0.5,ha='center',va='center',
		 bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.5})



#plt.legend(proxy,names,loc='upper right',prop={'size':12},frameon=False,ncol=2)
#CS = plt.contour(om,ol,-2.0*np.log(z),levels=[2.31, 6.17, 11.8] )
CS = plt.contour(om,ol,-2.0*np.log(z),levels=[2.31, 6.17],colors=['black','black'] )

plt.text(0.4,1.5,"R98 Discovery Sample", color='black',rotation=39,ha='center',va='center')
plt.text(0.33,0.75,"Pantheon", color='k',rotation=47,alpha=0.5,ha='center',va='center')#,
		 #bbox={'edgecolor':'r','facecolor':'r','alpha':0.5,'boxstyle':'round'})
plt.text(1.1,1.45,"RAISIN NIR (stat)", color='k',rotation=0,alpha=0.5,ha='center',va='center',
		 bbox={'edgecolor':'C1','facecolor':'C1','alpha':0.8,'boxstyle':'round'})

#plt.text(0.21,1.25,"Pantheon (Stat)", color='gray',rotation=40,alpha=0.5)


#CS2 = plt.contour(om,ol,z)

plt.figure(figsize=(20,10))
str2=colors
#colors=['black','green','blue','red','purple']
#'CMB', 'CMB+BAO','SN+CMB','SN'
#plt.savefig('plot_omol.png') 
g.export('plot_omol.png')
