import pylab as plt
plt.ion()
import numpy as np
import cosmo
from txtobj import txtobj
from astropy.stats import sigma_clipped_stats
import copy
#https://jiffyclub.github.io/palettable/

#maxmodelfile_raisin = 'nir_maxmodel_fitresults.txt'
#maxmodelfile_lowz = 'nir_maxmodel_fitresults_lowz.txt'
#maxmodelfile_stretchprior = ''

snanafile_raisin1 = 'output/fit_nir/PS1_RAISIN.FITRES.TEXT'
snanafile_raisin2 = 'output/fit_nir/DES_RAISIN.FITRES.TEXT'
snanafile_lowz_csp = 'output/fit_nir/CSP_RAISIN.FITRES.TEXT'
goodcids_raisin1 = np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',dtype=str,unpack=True)
goodcids_raisin2 = np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',dtype=str,unpack=True)
goodcids_lowz_csp = np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',dtype=str,unpack=True)

distkeys = ['Ymu','Jmu','Hmu']
disterrkeys = ['Ymuerr','Jmuerr','Hmuerr']
        
def hubbleplot():
    plt.rcParams['figure.figsize'] = (8,4)
    plt.rcParams['font.size'] = 13

    plt.clf()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
    
    ax1p1 = plt.axes([0.1,0.45,0.33,0.45])
    ax1p2 = plt.axes([0.48,0.45,0.33,0.45],sharey=ax1p1)
    ax2p1 = plt.axes([0.1,0.13,0.33,0.32])
    ax2p2 = plt.axes([0.48,0.13,0.33,0.32],sharey=ax2p1)
    ax3 = plt.axes([0.81,0.13,0.1,0.32],sharey=ax2p1)
    ax3.yaxis.tick_right()
    ax3.set_xticks([])

    import palettable

    from palettable.colorbrewer.qualitative import Dark2_8 as palettable_color
    for ax in [ax1p1,ax1p2,ax2p1,ax2p2,ax3]:
        ax.set_prop_cycle('color', palettable_color.mpl_colors)
        
    ax1p1.spines['right'].set_visible(False)
    ax1p2.spines['left'].set_visible(False)
    ax2p1.spines['right'].set_visible(False)
    ax2p2.spines['left'].set_visible(False)
    
    # get the distance moduli
    fr = txtobj(snanafile_raisin1,fitresheader=True)
    fr2 = txtobj(snanafile_raisin2,fitresheader=True)
    frlowz = txtobj(snanafile_lowz_csp,fitresheader=True)
    iGood = np.array([],dtype=int)
    for j,i in enumerate(fr.CID):
        if i in goodcids_raisin1: iGood = np.append(iGood,j)
    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iGood]
    iGood = np.array([],dtype=int)
    for j,i in enumerate(fr2.CID):
        if i in goodcids_raisin2: iGood = np.append(iGood,j)
    for k in fr2.__dict__.keys():
        fr2.__dict__[k] = fr2.__dict__[k][iGood]
    iGood = np.array([],dtype=int)
    for j,i in enumerate(frlowz.CID):
        if i in goodcids_lowz_csp: iGood = np.append(iGood,j)
    for k in frlowz.__dict__.keys():
        frlowz.__dict__[k] = frlowz.__dict__[k][iGood]

        
    for k in fr.__dict__.keys():
        fr.__dict__[k] = np.append(fr.__dict__[k],fr2.__dict__[k])

        
    print('%i RAISIN SNe'%len(fr.zHD))
    fr = getsnsurvey(fr)

    frlowz = getsnsurvey(frlowz)

    fr.mures = fr.DLMAG - cosmo.mu(fr.zHD)
    frlowz.mures = frlowz.DLMAG - cosmo.mu(frlowz.zHD)
    avgresid = sigma_clipped_stats(np.append(fr.DLMAG-cosmo.mu(fr.zHD),
                                             frlowz.DLMAG-cosmo.mu(frlowz.zHD)))[1]
    fr.DLMAG -= avgresid
    fr.mures -= avgresid
    frlowz.DLMAG -= avgresid
    frlowz.mures -= avgresid
    
    # plot the hubble resids
    iR1 = fr.IDSURVEY == 15
    iR2 = fr.IDSURVEY == 10
    ax1p1.errorbar(frlowz.zHD,frlowz.DLMAG,yerr=frlowz.DLMAGERR,fmt='o',
                   label='low-$z$',color=palettable_color.mpl_colors[0])
    ax1p1.errorbar(fr.zHD[iR1],fr.DLMAG[iR1],yerr=fr.DLMAGERR[iR1],fmt='o',
                   label='RAISIN1',color=palettable_color.mpl_colors[1])
    ax1p1.errorbar(fr.zHD[iR2],fr.DLMAG[iR2],yerr=fr.DLMAGERR[iR2],fmt='o',
                   label='RAISIN2',color=palettable_color.mpl_colors[2])

    ax1p2.errorbar(fr.zHD[iR1],fr.DLMAG[iR1],yerr=fr.DLMAGERR[iR1],fmt='o',
                   label='RAISIN1',color=palettable_color.mpl_colors[1])
    ax1p2.errorbar(fr.zHD[iR2],fr.DLMAG[iR2],yerr=fr.DLMAGERR[iR2],fmt='o',
                   label='RAISIN2',color=palettable_color.mpl_colors[2])
    plotlcdm(ax1p1)
    plotlcdm(ax1p2)
    ax2p1.errorbar(frlowz.zHD,frlowz.mures,yerr=frlowz.DLMAGERR,fmt='o',
                   label='low-$z$',color=palettable_color.mpl_colors[0])
    ax2p2.errorbar(fr.zHD[iR1],fr.mures[iR1],yerr=fr.DLMAGERR[iR1],fmt='o',
                   color=palettable_color.mpl_colors[1])
    ax2p2.errorbar(fr.zHD[iR2],fr.mures[iR2],yerr=fr.DLMAGERR[iR2],fmt='o',
                   color=palettable_color.mpl_colors[2])
    ax2p1.axhline(0,color='k',lw=2)
    ax2p2.axhline(0,color='k',lw=2)
    
    # making the plots look good
    ax1p1.set_xlim([0.01,0.05])
    ax1p2.set_xlim([0.2,0.7])
    ax2p1.set_xlim([0.01,0.05])
    ax2p2.set_xlim([0.2,0.7])
    ax1p1.legend(loc='upper left')

    ax3.hist(np.append(fr.mures,frlowz.mures[frlowz.zHD > 0.01]),bins=np.arange(-1,1,0.1),
             orientation='horizontal',color='r',density=True,alpha=0.7,ec='r',lw=2,label='NIR')
    frpan = txtobj('hlsp_ps1cosmo_panstarrs_gpc1_all_model_v1_ancillary-g10.fitres.txt',fitresheader=True)
    ax3.hist(frpan.MURES,bins=np.arange(-1,1,0.1),
             orientation='horizontal',color='b',density=True,alpha=1.0,ec='k',lw=2,label='Optical\n(Pantheon)',histtype='step')
    ax3.legend(prop={'size':10},bbox_to_anchor=(0.0,1.6),bbox_transform=ax3.transAxes,loc='upper left')

    from matplotlib.ticker import NullFormatter
    for ax in [ax1p1,ax1p2]:
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())

        ax.set_xticks([])
        ax.set_xlabel(r'$z_{\mathrm{CMB}}$',fontsize=15,labelpad=0)

        ax.set_ylabel('$\mu$ (mag)',fontsize=15)
        ax.set_ylim([33,45])
        ax.tick_params(top="off",bottom="on",direction="inout",length=8, width=2)
    for ax in [ax2p1,ax2p2]:
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())

        ax.set_xlabel(r'$z_{\mathrm{CMB}}$',fontsize=15,labelpad=0)

        ax.set_ylabel('$\mu - \mu_{\Lambda CDM}$ (mag)',fontsize=15,labelpad=0)
        ax.set_ylim([-1,1])
        ax.set_ylim([-0.5,0.5])
    
    ax.tick_params(top="off",bottom="on",direction="inout",length=8, width=2)

        
    ax1p2.set_ylabel('')
    ax2p2.set_ylabel('')
    ax1p1.set_xticks([0.01,0.03])
    ax1p1.set_xticklabels(['0.01','0.03'])
    ax1p2.set_xticks([0.3,0.5,0.7])
    ax1p2.set_xticklabels(['0.3','0.5','0.7'])
    ax2p1.set_xticks([0.01,0.03])
    ax2p1.set_xticklabels(['0.01','0.03'])
    ax2p2.set_xticks([0.3,0.5,0.7])
    ax2p2.set_xticklabels(['0.3','0.5','0.7'])
    ax1p2.set_xlim([0.2,0.8])
    ax2p2.set_xlim([0.2,0.8])
    ax1p1.set_yticks([35,37.5,40,42.5,45])
    ax1p1.set_xticklabels([])
    
    ax1p2.yaxis.tick_right()
    plt.setp(ax1p2.get_yticklabels(), visible=False)
    
    ax2p2.yaxis.tick_right()
    plt.setp(ax2p2.get_yticklabels(), visible=False)


    
    d = .015
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1p1.transAxes, color='k', clip_on=False)
    ax1p1.plot((1-d,1+d), (-d,+d), **kwargs)
    ax1p1.plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=ax1p2.transAxes)
    ax1p2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax1p2.plot((-d,+d), (-d,+d), **kwargs)

    kwargs = dict(transform=ax2p1.transAxes, color='k', clip_on=False)
    ax2p1.plot((1-d,1+d), (-d,+d), **kwargs)
    ax2p1.plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=ax2p2.transAxes)  # switch to the bottom axes
    ax2p2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax2p2.plot((-d,+d), (-d,+d), **kwargs)

    
    import pdb; pdb.set_trace()
        
def getsnsurvey(fr):

    fr.survey = []
    for i,snid in enumerate(fr.CID):
        if snid.startswith('DES') or snid.startswith('SNA'):
            fr.survey += ['RAISIN2']
        else:
            fr.survey += ['RAISIN1']
    fr.survey = np.array(fr.survey)

    return fr
    
def hubbleplotpars(ax):

    ax.set_xscale('log')
    ax.set_xlabel(r'$z_{\mathrm{CMB}}$',fontsize=15,labelpad=0)
    ax.set_xlim([0.01,0.75])

    ax.set_ylabel('$\mu$ (mag)',fontsize=15)
    ax.set_ylim([33,45])
    ax.tick_params(top="off",bottom="on",direction="inout",length=8, width=2)

    ax.set_xticklabels([])

def hubbleresidplotpars(ax):

    ax.set_xscale('log')
    ax.set_xlabel(r'$z_{\mathrm{CMB}}$',fontsize=15,labelpad=0)
    ax.set_xlim([0.01,0.75])

    ax.set_ylabel('$\mu - \mu_{\Lambda CDM}$ (mag)',fontsize=15,labelpad=0)
    ax.set_ylim([-1,1])
    ax.set_ylim([-0.5,0.5])

    ax.tick_params(top="off",bottom="on",direction="inout",length=8, width=2)

    
def plotlcdm(ax):

    zrange = np.arange(0,1,0.01)
    ax.plot(zrange,cosmo.mu(zrange),color='k',lw=2)

def weighted_avg_and_err(values, weights):
    """
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, np.sqrt(variance/len(values)))         

if __name__ == "__main__":
    hubbleplot()
