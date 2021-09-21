import getdist.plots as gplot
import matplotlib.pyplot as plt
import numpy as np
import plotsetup
import scipy
import scipy.signal

def density_contour_data(x, y, covariance_factor=None, n_bins=None, n_sigma=(1, 2)):
    r"""Generate the data for a plot with confidence contours of the density
    of points (useful for MCMC analyses).

    Parameters:

    - `x`, `y`: lists or numpy arrays with the x and y coordinates of the points
    - `covariance_factor`: optional, numerical factor to tweak the smoothness
    of the contours. If not specified, estimated using Scott's/Silverman's rule.
    The factor should be between 0 and 1; larger values means more smoothing is
    applied.
    - n_bins: number of bins in the histogram created as an intermediate step.
      this usually does not have to be changed.
    - n_sigma: integer or iterable of integers specifying the contours
      corresponding to the number of sigmas to be drawn. For instance, the
      default (1, 2) draws the contours containing approximately 68 and 95%
      of the points, respectively.
    """
    if n_bins is None:
        n_bins = min(10*int(np.sqrt(len(x))), 200)
    f_binned, x_edges, y_edges = np.histogram2d(x, y, normed=True, bins=n_bins)
    x_centers = (x_edges[:-1] + x_edges[1:])/2.
    y_centers = (y_edges[:-1] + y_edges[1:])/2.
    x_mean = np.mean(x_centers)
    y_mean = np.mean(y_centers)
    dataset = np.vstack([x, y])

    d = 2 # no. of dimensions

    if covariance_factor is None:
        # Scott's/Silverman's rule
        n = len(x) # no. of data points
        _covariance_factor = n**(-1/6.)
    else:
        _covariance_factor = covariance_factor

    cov = np.cov(dataset) * _covariance_factor**2
    gaussian_kernel = scipy.stats.multivariate_normal(mean=[x_mean, y_mean], cov=cov)

    x_grid, y_grid = np.meshgrid(x_centers, y_centers)
    xy_grid = np.vstack([x_grid.ravel(), y_grid.ravel()])
    f_gauss = gaussian_kernel.pdf(xy_grid.T)
    f_gauss = np.reshape(f_gauss, (len(x_centers), len(y_centers))).T

    f = scipy.signal.fftconvolve(f_binned, f_gauss, mode='same').T
    f = f/f.sum()

    def find_confidence_interval(x, pdf, confidence_level):
        return pdf[pdf > x].sum() - confidence_level
    def get_level(n):
        return scipy.optimize.brentq(find_confidence_interval, 0., 1.,
                                     args=(f.T, confidence_level(n)))
    #if isinstance(n_sigma, Number):
	#levels = [get_level(n_sigma)]
    #else:
    levels = [get_level(m) for m in sorted(n_sigma)]

    # replace negative or zero values by a tiny number before taking the log
    f[f <= 0] = 1e-32
    # convert probability to -2*log(probability), i.e. a chi^2
    f = -2*np.log(f)
    # convert levels to chi^2 and make the mode equal chi^2=0
    levels = list(-2*np.log(levels) - np.min(f))
    f = f - np.min(f)

    return {'x': x_grid, 'y': y_grid, 'z': f, 'levels': levels}

def area(vs):
    a = 0
    x0,y0 = vs[0]
    for [x1,y1] in vs[1:]:
        dx = x1-x0
        dy = y1-y0
        a += 0.5*(y0*dx - x0*dy)
        x0 = x1
        y0 = y1
    return a


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
g = gplot.getPlotter(chain_dir = '/project2/rkessler/SURVEYS/PS1MD/USERS/djones/RAISIN/cosmo')
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

roots = ['chains_raisin/noSN_omol', 'chains_raisin/DS17_ALL_omol_alone','chains_raisin/DS17_ALL_omol']
#roots = ['chains2/noSN_omol', 'chains2/DS17_ALL_omol_alone','chains2/DS17_ALL_omol']
#roots = ['chains2/noSN_omol','chains2/noSN_omol','chains2/noSN_omol','chains2/noSN_omol']
names=['CMB','SN','SN+CMB']
#names=['SN','SN(stat)','SN+CMB']

roots=['chains_raisin/RAISIN_all_ocdm','chains_raisin/DS17_ALL_omol_alone']
colors=['#C1292E','RAISIN (stat+sys)']
#colors=['teal','k','blue','purple']
#'#778899'
#colors=['#235789','#C1292E','#778899','#020100','#F1D302','#FDFFFC']
colors=['C1','red']
x=np.arange(0,2,.1)

#import pdb; pdb.set_trace()
#sagd = gplot.SampleAnalysisGetDist(g.plot_data)
#samples = g.sampleAnalyser.samplesForRoot(roots[0]).getParams()
#density = sagd.load_2d(samples,'omegam','omegal*')
#density.contours

g.plot_2d(roots,'omegam','omegal*',
          filled=[True,True,False,True],colors=colors)

line, = plt.plot(x, 1-x, lw=2,color='red',linestyle='--',alpha=1.0)
line, = plt.plot(x, 0.5*x, lw=2,color='blue',linestyle='--',alpha=1.0)
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

plt.text(1.2,0.67,"Accelerating Universe", color='blue',rotation=19,alpha=1.0,ha='center',va='center')
plt.text(1.22,0.53,"Decelerating Universe", color='blue',rotation=19,alpha=1.0,ha='center',va='center')

plt.text(0.65,0.24,"Flat Universe", color='red',rotation=325,alpha=1.0,ha='center',va='center',
		 bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.5})



#plt.legend(proxy,names,loc='upper right',prop={'size':12},frameon=False,ncol=2)
#CS = plt.contour(om,ol,-2.0*np.log(z),levels=[2.31, 6.17, 11.8] )
CS = plt.contour(om,ol,-2.0*np.log(z),levels=[2.31, 6.17],colors=['black','black'] )
levels=[2.31,6.17]
for i in range(len([2.31, 6.17])):
    contour = CS.collections[i]
    vs = contour.get_paths()[0].vertices
    # Compute area enclosed by vertices.
    a = area(vs)
    print("r = " + str(levels[i]) + ": a =" + str(a))

#ombins=np.arange(0,2.6,0.01)
#olbins=np.arange(-1,3.0,0.01)
#samples = g.sampleAnalyser.samplesForRoot(roots[0]).getParams()
#import flavio
#from flavio.statistics.functions import delta_chi2, confidence_level
#from numbers import Number
#outdict = density_contour_data(samples.omegam,samples.omegal)
#CS = plt.contour(outdict['x'],outdict['y'],outdict['z'],levels=levels,colors=['orange','black'])

#for i in range(len(levels)):
#    contour = CS.collections[i]
#    vs = contour.get_paths()[0].vertices
#    # Compute area enclosed by vertices.
#    a = area(vs)
#    print("r = " + str(levels[i]) + ": a =" + str(a))

#import pdb; pdb.set_trace()
#CS = plt.contour(xedges,yedges,H,levels=[0.68,0.95],colors=['orange','orange'])

plt.text(0.4,1.5,"R98 Discovery Sample", color='black',rotation=39,ha='center',va='center')
plt.text(0.33,0.75,"Pantheon", color='k',rotation=47,alpha=1.0,ha='center',va='center',fontweight='bold')#,
		 #bbox={'edgecolor':'r','facecolor':'1.0','alpha':0.3,'boxstyle':'round'})#,
		 #bbox={'edgecolor':'r','facecolor':'r','alpha':0.5,'boxstyle':'round'})
plt.text(0.7,1.05,"RAISIN NIR", color='k',rotation=0,alpha=1.0,ha='center',va='center',fontweight='bold',
		 bbox={'edgecolor':'1.0','facecolor':'1.0','alpha':0.3,'boxstyle':'round'})

#plt.text(0.21,1.25,"Pantheon (Stat)", color='gray',rotation=40,alpha=0.5)


#CS2 = plt.contour(om,ol,z)

plt.figure(figsize=(20,10))
str2=colors
#colors=['black','green','blue','red','purple']
#'CMB', 'CMB+BAO','SN+CMB','SN'
#plt.savefig('plot_omol.png') 
g.export('plot_omol.png')
