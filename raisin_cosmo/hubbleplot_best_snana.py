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
snanafile_lowz = 'output/fit_nir/LOWZ_RAISIN.FITRES.TEXT'
snanafile_lowz_csp = 'output/fit_nir/CSP_RAISIN.FITRES.TEXT'
snanafile_lowz_cfa = 'output/fit_nir/CfA_RAISIN.FITRES.TEXT'
maxmodelfile_raisin = '../snoopy_nir_pipeline/nir_maxmodel_fitresults.txt'



distkeys = ['Ymu','Jmu','Hmu']
disterrkeys = ['Ymuerr','Jmuerr','Hmuerr']

def hubbletable():
    fr = txtobj(maxmodelfile_raisin)
    fr.zcmb = fr.z[:]
    fr = getavgdist(fr)
    fr = getsnsurvey(fr)
    fr2 = mkcuts(copy.deepcopy(fr))

    idx = np.argsort(fr.z)
    for j,i in enumerate(fr.SNID[idx]):
        if fr.avgerr[idx][j] > 0.19: cutstr = 'Bad Host Sub.'
        elif i == 'PScF520107': cutstr = 'possible non-Ia'
        else: cutstr = '\\nodata'
        
        outline = "%s&%.3f&%.3f $\\pm$ %.3f&%.3f&%.3f&%s\\\\"%(
            i,fr.muavg[idx][j],fr.muavgerr[idx][j],fr.pkmjderr[idx][j],fr.avgerr[idx][j],cutstr)
        print(outline)
        
def add_optical_raisin_info(fr):
    
    opt_fr = txtobj('../snoopy_nir_pipeline/optical_fitresults.txt')
    fr_avgerr = txtobj(maxmodelfile_raisin)
    fr.opt_stretch,fr.opt_ebv,fr.opt_fitprob,fr.avgerr = \
        np.array([-99.0]*len(fr.CID)),np.array([-99.0]*len(fr.CID)),np.array([-99.0]*len(fr.CID)),np.array([-99.0]*len(fr.CID))

    for j,i in enumerate(fr.CID):
        if i == 'SNABELL370': fr.opt_stretch[j] = 1; fr.opt_ebv[j] = 0
        else:
            if not i in opt_fr.SNID: raise RuntimeError('missing optical fit for %s'%fr.CID[j])
            fr.opt_stretch[j] = opt_fr.st[opt_fr.SNID == i][0]
            fr.opt_ebv[j] = opt_fr.EBV[opt_fr.SNID == i][0]
            #fr.opt_fitprob[j] = opt_fr.FITPROB[opt_fr.CID == i][0]
            fr.avgerr[j] = fr_avgerr.avgerr[fr_avgerr.SNID == i][0]
        
    iBad = (fr.opt_ebv > 0.3) | (fr.opt_stretch < 0.8) | (fr.opt_stretch > 1.3) #| (fr.opt_fitprob < 0.001)
    print(fr.CID[iBad],fr.opt_stretch[iBad],fr.opt_ebv[iBad])
    
    iGood = (fr.opt_ebv < 0.3) & (fr.opt_stretch > 0.8) & (fr.opt_stretch < 1.3) #& (fr.opt_fitprob > 0.001)
    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iGood]

        
    return fr
        
def add_optical_lowz_info(fr):
    opt_fr = txtobj('output/fit_nir/CSP_RAISIN_optical.FITRES.TEXT',fitresheader=True)
    #opt_fr = txtobj('../snoopy_nir_pipeline/optical_lowz_fitresults.txt')
    fr.opt_stretch,fr.opt_ebv = np.array([-99.0]*len(fr.CID)),np.array([-99.0]*len(fr.CID))
    for j,i in enumerate(fr.CID):
        if not i in opt_fr.CID:
            print('missing optical fit for %s'%fr.CID[j])
            fr.opt_stretch[j] = 99
            fr.opt_ebv[j] = 99
            continue
            #raise RuntimeError('missing optical fit for %s'%fr.CID[j])
        
        fr.opt_stretch[j] = opt_fr.STRETCH[opt_fr.CID == i][0]
        fr.opt_ebv[j] = opt_fr.AV[opt_fr.CID == i][0]/3.1

    iGood = (fr.opt_ebv < 0.3) & (fr.opt_stretch > 0.8) & (fr.opt_stretch < 1.3)
    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iGood]
    return fr

def hubbletable_new():
    cutsdict = {'DES15C1nhv':'Chauvenet',
                'DES16E2cxw':'bad host sub.',
                'DES16E2rd':'bad host sub.',
                'DES16X1cpf':'bad host sub.',
                'PScA470041':'bad host sub.',
                'PScF520107':'possible non-Ia',
                'PScK450082':'bad host sub.',
                #'PScD500301':'Chauvenet',
                'PScD500100':'low $s_{BV}$'}
                #'SNABELL370 lazy}
    #snanafile_raisin1 = 'output/fit_nir/PS1_RAISIN.FITRES.TEXT'
    #snanafile_raisin2 = 'output/fit_nir/DES_RAISIN.FITRES.TEXT'
    frp = txtobj('output/fit_nir_sys/PS1_RAISIN/PS1_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    frd = txtobj('output/fit_nir_sys/DES_RAISIN/DES_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    lp = txtobj('output/fit_nir_sys/PS1_RAISIN/PS1_RAISIN/FITOPT000.LCPLOT.gz',fitresheader=True,rowprfx='OBS')
    ld = txtobj('output/fit_nir_sys/DES_RAISIN/DES_RAISIN/FITOPT000.LCPLOT.gz',fitresheader=True,rowprfx='OBS')

    frm = txtobj(maxmodelfile_raisin)
    fr = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT000_new.FITRES',
                fitresheader=True)
    idx = np.argsort(frm.z)
    for j,i in enumerate(frm.SNID[idx]):
        if i in cutsdict.keys():
            if i.startswith('PSc'):i2 = f'PS1-{i[-6:]}'
            else: i2 = i[:]
            outline = f"%s&{frm.z[idx][j]:.3f}&\\nodata&\\nodata&\\nodata&%.3f&%.3f&%s\\\\"%(
                i2,frm.pkmjderr[idx][j],frm.avgerr[idx][j],cutsdict[i])
            print(outline)
        elif i.startswith('SN'):
            pass
            #outline = "%s&\\nodata&%.3f&%.3f&\\nodata\\\\"%(
            #    i,frm.muavg[j],frm.muavgerr[j],frm.pkmjderr[j],frm.avgerr[j])
        elif i.startswith('PScC490521'):
            biascor = fr.DLMAG_biascor[fr.CID == i][0] #- frp.DLMAG[frp.CID == i][0]
            i2 = i.replace('PScC','PS1-')
            rest_bands = ','.join(np.unique(lp.BAND_REST[(lp.CID == i) & (lp.BAND_REST != '!')])[::-1])
            print(f"{i2}&{frm.z[idx][j]:.3f}&{frp.DLMAG[frp.CID == i][0]:.3f} $\\pm$ {frp.DLMAGERR[frp.CID == i][0]:.3f}&${rest_bands}$&{biascor:.3f}&{frm.pkmjderr[frm.SNID == i][0]:.3f}&{frm.avgerr[frm.SNID == i][0]:.3f}&\\nodata\\\\")       
        else:
            if i in frp.CID:
                try:
                    biascor = fr.DLMAG_biascor[fr.CID == i][0]# - frp.DLMAG[frp.CID == i][0]
                    i2 = f'PS1-{i[-6:]}'
                    rest_bands = ','.join(np.unique(lp.BAND_REST[(lp.CID == i) & (lp.BAND_REST != '!')])[::-1])
                    print(f"{i2}&{frm.z[idx][j]:.3f}&{frp.DLMAG[frp.CID == i][0]:.3f} $\\pm$ {frp.DLMAGERR[frp.CID == i][0]:.3f}&${rest_bands}$&{biascor:.3f}&{frm.pkmjderr[frm.SNID == i][0]:.3f}&{frm.avgerr[frm.SNID == i][0]:.3f}&\\nodata\\\\")
                except: import pdb; pdb.set_trace()
            elif i in frd.CID:
                biascor = fr.DLMAG_biascor[fr.CID == i][0] #- frd.DLMAG[frd.CID == i][0]
                rest_bands = ','.join(np.unique(ld.BAND_REST[(ld.CID == i) & (ld.BAND_REST != '!')])[::-1])
                print(f"{i}&{frm.z[idx][j]:.3f}&{frd.DLMAG[frd.CID == i][0]:.3f} $\\pm$ {frd.DLMAGERR[frd.CID == i][0]:.3f}&${rest_bands}$&{biascor:.3f}&{frm.pkmjderr[frm.SNID == i][0]:.3f}&{frm.avgerr[frm.SNID == i][0]:.3f}&\\nodata\\\\")

                
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
    #from palettable.wesanderson import Aquatic1_5 as palettable_color

    from palettable.colorbrewer.qualitative import Dark2_8 as palettable_color
    for ax in [ax1p1,ax1p2,ax2p1,ax2p2,ax3]:
        ax.set_prop_cycle('color', palettable_color.mpl_colors)

    #palettable_color = ['#66c2a5',
    #                   '#fc8d62',
    #   '#8da0cb']
        
    ax1p1.spines['right'].set_visible(False)
    ax1p2.spines['left'].set_visible(False)
    ax2p1.spines['right'].set_visible(False)
    ax2p2.spines['left'].set_visible(False)
    
    # get the distance moduli
    fr = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT000.FITRES',fitresheader=True)
    #fr = txtobj(snanafile_raisin1,fitresheader=True)
    #fr2 = txtobj(snanafile_raisin2,fitresheader=True)
    #for k in fr.__dict__.keys():
    #    fr.__dict__[k] = np.append(fr.__dict__[k],fr2.__dict__[k])
    #fr_avgerr = txtobj(maxmodelfile_raisin)
    #fr.avgerr = np.array([-99.0]*len(fr.CID))
    #for j,i in enumerate(fr.CID):
    #    fr.avgerr[j] = fr_avgerr.avgerr[fr_avgerr.SNID == i][0]

        
    print('%i RAISIN SNe before cuts'%len(fr.zCMB))
    #fr = getavgdist(fr)
    #fr = getsnsurvey(fr)
    #fr = mkcuts(fr)
    print('%i RAISIN SNe after cuts'%len(fr.zCMB))
    #fr = add_optical_raisin_info(fr)
    print('%i RAISIN SNe after optical cuts'%len(fr.zCMB))

    frlowz = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT000.FITRES',fitresheader=True)
    iCSP = frlowz.IDSURVEY == 5
    for k in frlowz.__dict__.keys():
       frlowz.__dict__[k] = frlowz.__dict__[k][iCSP]

    #frlowz = txtobj(snanafile_lowz_csp,fitresheader=True)
    #frlowz2 = txtobj(snanafile_lowz_cfa,fitresheader=True)
    #for k in frlowz.__dict__.keys():
    #   frlowz.__dict__[k] = np.append(frlowz.__dict__[k],frlowz2.__dict__[k])
    #import pdb; pdb.set_trace()
    
    print('%i low-z SNe before cuts'%len(frlowz.zCMB))
    frlowz = getlowzpars(frlowz)
    #frlowz = getavgdist(frlowz)
    frlowz = getsnsurvey(frlowz)
    frlowz = mkcuts(frlowz)
    #print('%i low-z SNe after cuts'%len(frlowz.zCMB))
    #frlowz = add_optical_lowz_info(frlowz)
    #print('%i low-z SNe after optical cuts'%len(frlowz.zCMB))

    fr.mures = fr.DLMAG - cosmo.mu(fr.zCMB)
    frlowz.mures = frlowz.DLMAG - cosmo.mu(frlowz.zcmb)
    avgresid = sigma_clipped_stats(np.append(fr.DLMAG-cosmo.mu(fr.zCMB),
                                             frlowz.DLMAG-cosmo.mu(frlowz.zcmb)))[1]
    fr.DLMAG -= avgresid
    fr.mures -= avgresid
    frlowz.DLMAG -= avgresid
    frlowz.mures -= avgresid
    
    # plot the hubble resids
    iR1 = fr.IDSURVEY == 15
    iR2 = fr.IDSURVEY == 10
    ax1p1.errorbar(frlowz.zcmb,frlowz.DLMAG,yerr=frlowz.DLMAGERR,fmt='o',
                   label='low-$z$',color=palettable_color.mpl_colors[0])
    ax1p1.errorbar(fr.zCMB[iR1],fr.DLMAG[iR1],yerr=fr.DLMAGERR[iR1],fmt='o',
                   label='RAISIN1',color=palettable_color.mpl_colors[1])
    ax1p1.errorbar(fr.zCMB[iR2],fr.DLMAG[iR2],yerr=fr.DLMAGERR[iR2],fmt='o',
                   label='RAISIN2',color=palettable_color.mpl_colors[2])

    ax1p2.errorbar(fr.zCMB[iR1],fr.DLMAG[iR1],yerr=fr.DLMAGERR[iR1],fmt='o',
                   label='RAISIN1',color=palettable_color.mpl_colors[1])
    ax1p2.errorbar(fr.zCMB[iR2],fr.DLMAG[iR2],yerr=fr.DLMAGERR[iR2],fmt='o',
                   label='RAISIN2',color=palettable_color.mpl_colors[2])
    plotlcdm(ax1p1)
    plotlcdm(ax1p2)
    ax2p1.errorbar(frlowz.zcmb,frlowz.mures,yerr=frlowz.DLMAGERR,fmt='o',
                   label='low-$z$',color=palettable_color.mpl_colors[0])
    ax2p2.errorbar(fr.zCMB[iR1],fr.mures[iR1],yerr=fr.DLMAGERR[iR1],fmt='o',
                   color=palettable_color.mpl_colors[1])
    ax2p2.errorbar(fr.zCMB[iR2],fr.mures[iR2],yerr=fr.DLMAGERR[iR2],fmt='o',
                   color=palettable_color.mpl_colors[2])
    ax2p1.axhline(0,color='k',lw=2)
    ax2p2.axhline(0,color='k',lw=2)
    
    # making the plots look good
    #hubbleplotpars(ax1p1)
    #hubbleplotpars(ax1p2)
    ax1p1.set_xlim([0.01,0.05])
    ax1p2.set_xlim([0.2,0.7])
    #hubbleresidplotpars(ax2p1)
    #hubbleresidplotpars(ax2p2)
    ax2p1.set_xlim([0.01,0.05])
    ax2p2.set_xlim([0.2,0.7])
    ax1p1.legend(loc='upper left')

    ax3.hist(np.append(fr.mures,frlowz.mures[frlowz.zcmb > 0.01]),bins=np.arange(-1,1,0.1),
             orientation='horizontal',color='r',density=True,alpha=0.7,ec='r',lw=2,label='NIR')#,histtype='step')# (low-$z$ + RAISIN)')
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

        
    #ax1p2.yaxis.set_ticks([])
    #ax2p2.yaxis.set_ticks([])
    ax1p2.set_ylabel('')
    ax2p2.set_ylabel('')
    ax1p1.set_xticks([0.01,0.03]) #,0.3,0.5,0.7])
    ax1p1.set_xticklabels(['0.01','0.03']) #,'0.3','0.5','0.7'])
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
    
    #ax1p1.yaxis.tick_left()
    #ax1p1.tick_params(labelright='off')
    ax1p2.yaxis.tick_right()
    plt.setp(ax1p2.get_yticklabels(), visible=False)
    
    #ax2p1.yaxis.tick_left()
    #ax2p1.tick_params(labelright='off')
    ax2p2.yaxis.tick_right()
    plt.setp(ax2p2.get_yticklabels(), visible=False)


    
    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1p1.transAxes, color='k', clip_on=False)
    ax1p1.plot((1-d,1+d), (-d,+d), **kwargs)
    ax1p1.plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=ax1p2.transAxes)  # switch to the bottom axes
    ax1p2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax1p2.plot((-d,+d), (-d,+d), **kwargs)

    kwargs = dict(transform=ax2p1.transAxes, color='k', clip_on=False)
    ax2p1.plot((1-d,1+d), (-d,+d), **kwargs)
    ax2p1.plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=ax2p2.transAxes)  # switch to the bottom axes
    ax2p2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax2p2.plot((-d,+d), (-d,+d), **kwargs)

    #print(fr.CID)
    
    import pdb; pdb.set_trace()
    
def mkcuts(fr):

    if 'ebv' in fr.__dict__.keys(): iCut = fr.ebv < 0.25
    elif 'AV' in fr.__dict__.keys(): iCut = fr.AV < 0.25*3.1
    if 'avgerr' in fr.__dict__.keys():
        iCut = (fr.avgerr < 0.19)
        print(fr.CID[fr.avgerr > 0.19])

    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iCut]

    iCut = (fr.CID != 'PScF520107') & (fr.CID != 'snf20080514-002') & (fr.CID != 'PScF510457') & (fr.CID != 'sn2005bo') # & (fr.CID != 'SNABELL370')
    iCut = (fr.CID != 'snf20080514-002') & (fr.CID != 'PScF510457') & (fr.CID != 'sn2005bo') # & (fr.CID != 'SNABELL370')
    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iCut]

    iCut = fr.DLMAGERR < 0.4
    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iCut]

    return fr
    
def getlowzpars(fr):
    lowzpars = txtobj('../snoopy_nir_pipeline/GPfit/lcparams_avelino19.txt')
    
    fr.zcmb = np.zeros(len(fr.zCMB))
    fr.ebv = np.zeros(len(fr.zCMB))
    for i in range(len(fr.CID)):
        try:
            fr.zcmb[i] = lowzpars.z_cmb[lowzpars.snid == 'SN'+fr.CID[i][2:]][0]
            fr.ebv[i] = lowzpars.EBVMW[lowzpars.snid == 'SN'+fr.CID[i][2:]][0]
        except:
            fr.zcmb[i] = fr.zCMB[i]
            fr.ebv[i] = fr.MWEBV[i]
            
    return fr
        
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
    #ax.set_xticks([0.01,0.1,0.3,0.5,0.7])
    #ax.set_xticklabels(['0.01','0.1','0.3','0.5','0.7'])
    ax.set_ylim([33,45])
    ax.tick_params(top="off",bottom="on",direction="inout",length=8, width=2)

    ax.set_xticklabels([])

def hubbleresidplotpars(ax):

    ax.set_xscale('log')
    ax.set_xlabel(r'$z_{\mathrm{CMB}}$',fontsize=15,labelpad=0)
    ax.set_xlim([0.01,0.75])

    ax.set_ylabel('$\mu - \mu_{\Lambda CDM}$ (mag)',fontsize=15,labelpad=0)
    ax.set_ylim([-1,1])
    #ax.set_xticks([0.01,0.1,0.3,0.5,0.7])
    #ax.set_xticklabels(['0.01','0.1','0.3','0.5','0.7'])
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
    #hubbletable_new()
    hubbleplot()
