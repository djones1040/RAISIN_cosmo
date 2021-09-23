#!/usr/bin/env python
# D. Jones - 1/27/21

from txtobj import txtobj
import numpy as np
import pylab as plt
plt.ion()
import cosmo
import scipy
import copy
import getmu

_nirfile = 'output/cosmo_fitres/RAISIN_combined_FITOPT000_nobiascor.FITRES'
_csprvfile = 'output/fit_optical/CSP_RAISIN_optnir_MWRV.FITRES.TEXT'
_ps1rvfile = 'output/fit_optical/PS1_RAISIN_optnir_MWRV.FITRES.TEXT'
_desrvfile = 'output/fit_optical/DES_RAISIN_optnir_MWRV.FITRES.TEXT'
_cspoptfile = 'output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT'
_ps1optfile = 'output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT'
_desoptfile = 'output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT'
#_csprvfile = 'output/fit_optical/CSP_RAISIN_SALT2.FITRES.TEXT'
#_ps1rvfile = 'output/fit_optical/PS1_RAISIN_SALT2.FITRES.TEXT'
#_desrvfile = 'output/fit_optical/DES_RAISIN_SALT2.FITRES.TEXT'

_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',dtype=str),
                            np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',dtype=str),
                            np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',dtype=str)))

# NIR vs. optical distances
# NIR vs. optical with RV floated

# let's pre-compute some useful crap

frnir = txtobj(_nirfile,fitresheader=True)
# here's normal opt+NIR
froptcsp = txtobj(_cspoptfile,fitresheader=True)
froptps1 = txtobj(_ps1optfile,fitresheader=True)
froptdes = txtobj(_desoptfile,fitresheader=True)
frrvcsp = txtobj(_csprvfile,fitresheader=True)
frrvps1 = txtobj(_ps1rvfile,fitresheader=True)
frrvdes = txtobj(_desrvfile,fitresheader=True)

idxcsp,idxps1,idxdes = np.array([],dtype=int),np.array([],dtype=int),np.array([],dtype=int)
for j,i in enumerate(frnir.CID):
    if i not in _goodcids: continue
    if i in froptcsp.CID:
        idxcsp = np.append(idxcsp,np.where(froptcsp.CID == i)[0][0])
    if i in froptps1.CID:
        idxps1 = np.append(idxps1,np.where(froptps1.CID == i)[0][0])
    if i in froptdes.CID:
        idxdes = np.append(idxdes,np.where(froptdes.CID == i)[0][0])
for k in froptcsp.__dict__.keys():
    froptcsp.__dict__[k] = np.concatenate(
        (froptcsp.__dict__[k][idxcsp],froptps1.__dict__[k][idxps1],
        froptdes.__dict__[k][idxdes]))
_fropt = froptcsp

# here's opt+NIR with low R_V
idxcsp,idxps1,idxdes = np.array([],dtype=int),np.array([],dtype=int),np.array([],dtype=int)
for j,i in enumerate(frnir.CID):
    if i not in _goodcids: continue
    if i in frrvcsp.CID:
        idxcsp = np.append(idxcsp,np.where(frrvcsp.CID == i)[0][0])
    if i in frrvps1.CID:
        idxps1 = np.append(idxps1,np.where(frrvps1.CID == i)[0][0])
    if i in frrvdes.CID:
        idxdes = np.append(idxdes,np.where(frrvdes.CID == i)[0][0])

if 'PKMJDINI' in frrvcsp.__dict__.keys():
    del frrvcsp.__dict__['PKMJDINI']
for k in frrvcsp.__dict__.keys():
    frrvcsp.__dict__[k] = np.concatenate(
        (frrvcsp.__dict__[k][idxcsp],frrvps1.__dict__[k][idxps1],
        frrvdes.__dict__[k][idxdes]))
_frrv = frrvcsp

# here's NIR with optical parameters
_frnirmod = txtobj('st_av_corr_mags.txt')

# and finally SALT2
frscsp = txtobj('output/fit_optical/CSP_RAISIN_SALT3.FITRES.TEXT',fitresheader=True)
frsps1 = txtobj('output/fit_optical/PS1_RAISIN_SALT3.FITRES.TEXT',fitresheader=True)
frsdes = txtobj('output/fit_optical/DES_RAISIN_SALT3.FITRES.TEXT',fitresheader=True)
idxcsp,idxps1,idxdes = np.array([],dtype=int),np.array([],dtype=int),np.array([],dtype=int)
for j,i in enumerate(frnir.CID):
    if i not in _goodcids: continue
    if i in frscsp.CID:
        idxcsp = np.append(idxcsp,np.where(frscsp.CID == i)[0][0])
    if i in frsps1.CID:
        idxps1 = np.append(idxps1,np.where(frsps1.CID == i)[0][0])
    if i in frsdes.CID:
        idxdes = np.append(idxdes,np.where(frsdes.CID == i)[0][0])

for k in frscsp.__dict__.keys():
    frscsp.__dict__[k] = np.concatenate(
        (frscsp.__dict__[k][idxcsp],frsps1.__dict__[k][idxps1],
        frsdes.__dict__[k][idxdes]))
_frs = frscsp
#import pdb; pdb.set_trace()

def main():
    plt.subplots_adjust(hspace=0)
    
    frnir = txtobj(_nirfile,fitresheader=True)
    froptcsp = txtobj(_cspoptfile,fitresheader=True)
    froptps1 = txtobj(_ps1optfile,fitresheader=True)
    froptdes = txtobj(_desoptfile,fitresheader=True)
    frrvcsp = txtobj(_csprvfile,fitresheader=True)
    frrvps1 = txtobj(_ps1rvfile,fitresheader=True)
    frrvdes = txtobj(_desrvfile,fitresheader=True)

    idxcsp,idxps1,idxdes = np.array([],dtype=int),np.array([],dtype=int),np.array([],dtype=int)
    for j,i in enumerate(frnir.CID):
        if i in froptcsp.CID:
            idxcsp = np.append(idxcsp,np.where(froptcsp.CID == i)[0][0])
        if i in froptps1.CID:
            idxps1 = np.append(idxps1,np.where(froptps1.CID == i)[0][0])
        if i in froptdes.CID:
            idxdes = np.append(idxdes,np.where(froptdes.CID == i)[0][0])
    for k in froptcsp.__dict__.keys():
        froptcsp.__dict__[k] = np.concatenate(
            (froptcsp.__dict__[k][idxcsp],froptps1.__dict__[k][idxps1],
            froptdes.__dict__[k][idxdes]))
    fropt = froptcsp

    idxcsp,idxps1,idxdes = np.array([],dtype=int),np.array([],dtype=int),np.array([],dtype=int)
    for j,i in enumerate(frnir.CID):
        if i in frrvcsp.CID:
            idxcsp = np.append(idxcsp,np.where(frrvcsp.CID == i)[0][0])
        if i in frrvps1.CID:
            idxps1 = np.append(idxps1,np.where(frrvps1.CID == i)[0][0])
        if i in frrvdes.CID:
            idxdes = np.append(idxdes,np.where(frrvdes.CID == i)[0][0])

    if 'PKMJDINI' in frrvcsp.__dict__.keys():
        del frrvcsp.__dict__['PKMJDINI']
    for k in frrvcsp.__dict__.keys():
        frrvcsp.__dict__[k] = np.concatenate(
            (frrvcsp.__dict__[k][idxcsp],frrvps1.__dict__[k][idxps1],
            frrvdes.__dict__[k][idxdes]))
    frrv = frrvcsp
    #import getmu
    #frrv = getmu.getmu(frrv)
    #frrv.DLMAG = frrv.mu
    #frrv.DLMAGERR = frrv.muerr

    frrv.DLMAG[frrv.HOST_LOGMASS > 10] += 0.04
    frrv.DLMAG[frrv.HOST_LOGMASS < 10] -= 0.04
    fropt.DLMAG[fropt.HOST_LOGMASS > 10] += 0.04
    fropt.DLMAG[fropt.HOST_LOGMASS < 10] -= 0.04
    
    ax1,ax2 = plt.subplot(211),plt.subplot(212)#,plt.subplot(313)
    delmuopt = np.median(fropt.DLMAG[fropt.zHD > 0.15]-cosmo.mu(fropt.zHD[fropt.zHD > 0.15])) - \
               np.median(fropt.DLMAG[fropt.zHD < 0.15]-cosmo.mu(fropt.zHD[fropt.zHD < 0.15]))
    delmurv = np.median(frrv.DLMAG[frrv.zHD > 0.15]-cosmo.mu(frrv.zHD[frrv.zHD > 0.15])) - \
               np.median(frrv.DLMAG[frrv.zHD < 0.15]-cosmo.mu(frrv.zHD[frrv.zHD < 0.15]))
    delmunir = np.median(frnir.DLMAG[frnir.zHD > 0.15]-cosmo.mu(frnir.zHD[frnir.zHD > 0.15])) - \
               np.median(frnir.DLMAG[frnir.zHD < 0.15]-cosmo.mu(frnir.zHD[frnir.zHD < 0.15]))
    #import pdb; pdb.set_trace()
    fropt.DLMAG -= np.median(fropt.DLMAG-cosmo.mu(fropt.zHD))
    frrv.DLMAG -= np.median(frrv.DLMAG-cosmo.mu(frrv.zHD))
    frnir.DLMAG -= np.median(frnir.DLMAG-cosmo.mu(frnir.zHD))    
    ax1.errorbar(fropt.zHD,fropt.DLMAG-cosmo.mu(fropt.zHD),yerr=fropt.DLMAGERR,fmt='*',color='b',label=f'Optical+NIR, $R_V = 3.1$, med. $\Delta\mu = {delmuopt:.3f}$')
    ax1.errorbar(frrv.zHD,frrv.DLMAG-cosmo.mu(frrv.zHD),yerr=frrv.DLMAGERR,fmt='D',color='lightblue',label=f'Optical+NIR, $R_V = 1.5$, med. $\Delta\mu = {delmurv:.3f}$')
    ax1.errorbar(frnir.zHD,frnir.DLMAG-cosmo.mu(frnir.zHD),yerr=frnir.DLMAGERR,fmt='o',color='r',label=f'NIR, med. $\Delta\mu = {delmunir:.3f}$')
    ax1.axhline(0,color='k',lw=2)
    ax1.legend(loc='upper left',prop={'size':7})
    ax1.set_ylim([-1.0,1.0])
    
    ax2.errorbar(fropt.zHD,fropt.AV,yerr=fropt.AVERR,fmt='o',color='b',label='$A_V$ $(R_V = 2.0)$')
    A = np.vander(fropt.zHD, 2)
    C = np.diag(fropt.AVERR * fropt.AVERR)
    ATA = np.dot(A.T, A / (fropt.AVERR ** 2)[:, None])
    cov = np.linalg.inv(ATA)
    w = np.linalg.solve(ATA, np.dot(A.T, fropt.AV / fropt.AVERR ** 2))
    print("Least-squares estimates:")
    print("m = {0:.3f} ± {1:.3f}".format(w[0], np.sqrt(cov[0, 0])))
    print("b = {0:.3f} ± {1:.3f}".format(w[1], np.sqrt(cov[1, 1])))
    
    def objective(x, a, b):
	    return a * x + b
    z = np.linspace(0.0,0.7,100)
    av_line = objective(z, w[0], w[1])
    av_upper = objective(z, w[0]+np.sqrt(cov[0, 0]), w[1])
    av_lower = objective(z, w[0]-np.sqrt(cov[0, 0]), w[1])
    ax2.plot(z,av_line,color='0.7',ls='--',label=f'slope = {w[0]:.3f} ± {np.sqrt(cov[0, 0]):.3f}')
    ax2.fill_between(z,av_lower,av_upper,color='0.7',alpha=0.5)
    ax2.axhline(0.0,color='k',lw=2)
    ax2.legend(loc='upper left',prop={'size':7})

    ax2.set_xlabel('$z_{CMB}$',fontsize=15)
    ax1.set_ylabel('$\mu - \mu_{\Lambda CDM}$',fontsize=15)
    ax1.xaxis.set_ticklabels([])
    ax2.set_ylabel('$A_V$ ($R_V = 2.0$)',fontsize=15)
    #ax2.set_ylabel('$R_V$',fontsize=15)

    for ax in [ax1,ax2]:
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        ax.set_xlim([0.0,0.7])
    
    import pdb; pdb.set_trace()
    return
    
    ax3.errorbar(frrv.zHD,frrv.RV,yerr=frrv.RVERR,fmt='o',color='b',label='$R_V$')
    A = np.vander(fropt.zHD, 2)
    C = np.diag(fropt.RVERR * fropt.RVERR)
    ATA = np.dot(A.T, A / (fropt.RVERR ** 2)[:, None])
    cov = np.linalg.inv(ATA)
    w = np.linalg.solve(ATA, np.dot(A.T, fropt.RV / fropt.RVERR ** 2))
    print("Least-squares estimates:")
    print("m = {0:.3f} ± {1:.3f}".format(w[0], np.sqrt(cov[0, 0])))
    print("b = {0:.3f} ± {1:.3f}".format(w[1], np.sqrt(cov[1, 1])))
    
    def objective(x, a, b):
	    return a * x + b
    z = np.linspace(0.0,0.7,100)
    rv_line = objective(z, w[0], w[1])
    rv_upper = objective(z, w[0]+np.sqrt(cov[0, 0]), w[1])
    rv_lower = objective(z, w[0]-np.sqrt(cov[0, 0]), w[1])
    ax3.plot(z,rv_line,color='0.7',ls='--',label=f'slope = {w[0]:.3f} ± {np.sqrt(cov[0, 0]):.3f}')
    ax3.fill_between(z,rv_lower,rv_upper,color='0.7',alpha=0.5)
    ax3.axhline(2.02,color='k',lw=2,label='CSP med. $R_V$')
    ax3.legend()

def newfig():
    plt.rcParams['figure.figsize'] = (8,8)
    from matplotlib.gridspec import GridSpec
    fig = plt.figure()
    gs = GridSpec(4, 4, figure=fig)
    
    plt.subplots_adjust(wspace=0.7,hspace=0)
    ax1 = fig.add_subplot(gs[0,0:3])
    ax2 = fig.add_subplot(gs[1,0:3])
    ax3 = fig.add_subplot(gs[2,0:3])
    ax4 = fig.add_subplot(gs[3,0:3])
    ax1hist = fig.add_subplot(gs[0,3])
    ax2hist = fig.add_subplot(gs[1,3])
    ax3hist = fig.add_subplot(gs[2,3])
    ax4hist = fig.add_subplot(gs[3,3])

    frbase = txtobj(_nirfile,fitresheader=True)
    frbase.resid = frbase.DLMAG - cosmo.mu(frbase.zHD)
    mass_step_approx = np.average(frbase.resid[frbase.HOST_LOGMASS > 10],weights=1/frbase.DLMAGERR[frbase.HOST_LOGMASS > 10]**2.)-\
        np.average(frbase.resid[frbase.HOST_LOGMASS < 10],weights=1/frbase.DLMAGERR[frbase.HOST_LOGMASS < 10]**2.)
    print(mass_step_approx)
    frbase.resid[frbase.HOST_LOGMASS > 10] += np.abs(mass_step_approx)/2.
    frbase.resid[frbase.HOST_LOGMASS < 10] -= np.abs(mass_step_approx)/2.

    
    frbase.resid -= np.median(frbase.resid)
    frs = getmu.mkcuts(copy.deepcopy(_frs))
    #frs = copy.deepcopy(_frs)
    
    for ax,axhist,frvar,label in zip(
            [ax1,ax2,ax3,ax4],
            [ax1hist,ax2hist,ax3hist,ax4hist],
            [_fropt,_frrv,_frs,_frnirmod],
            ['Optical+NIR ($R_V = 3.1$)',
             'Optical+NIR ($R_V = 2.0$)',
             'Optical with SALT2',
             'Optical+NIR $s_{BV}$, NIR dist.']):

        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        axhist.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        
        if 'SALT2' in label:
            frvar = getmu.getmu(frvar)
            frvar = getmu.mkcuts(frvar)
            frvar.resid = frvar.mures
            frvar.DLMAGERR = frvar.muerr
            
        frbase_matched = copy.deepcopy(frbase)
        frbase_matched.match_to_col('CID',frs.CID)
        frbase_matched.match_to_col('CID',frvar.CID)
        frvar.match_to_col('CID',frbase_matched.CID)

        if 'resid' not in frvar.__dict__.keys():
            frvar.resid = frvar.DLMAG - cosmo.mu(frvar.zHD)
            
        mass_step_approx = np.average(frvar.resid[frvar.HOST_LOGMASS > 10],weights=1/frvar.DLMAGERR[frvar.HOST_LOGMASS > 10]**2.)-\
            np.average(frvar.resid[frvar.HOST_LOGMASS < 10],weights=1/frvar.DLMAGERR[frvar.HOST_LOGMASS < 10]**2.)
        print(mass_step_approx)
            
        frvar.resid[frvar.HOST_LOGMASS > 10] += np.abs(mass_step_approx)/2.
        frvar.resid[frvar.HOST_LOGMASS < 10] -= np.abs(mass_step_approx)/2.
            
        frvar.resid -= np.median(frvar.resid)
        #if 'SALT2' in label: frvar.resid[frvar.zHD > 0.1] -= 0.073
        #else: frvar.resid[frvar.zHD > 0.1] -= 0.019

        diff_lowz,differr_lowz = weighted_avg_and_err(
            frvar.resid[frvar.zHD < 0.1]-frbase_matched.resid[frvar.zHD < 0.1],
            1/(frvar.DLMAGERR[frvar.zHD < 0.1]**2.+frbase_matched.DLMAGERR[frvar.zHD < 0.1]**2.))
        diff_highz,differr_highz = weighted_avg_and_err(
            frvar.resid[frvar.zHD > 0.1]-frbase_matched.resid[frvar.zHD > 0.1],
            1/(frvar.DLMAGERR[frvar.zHD > 0.1]**2.+frbase_matched.DLMAGERR[frvar.zHD > 0.1]**2.))
        avgdiff = diff_highz - diff_lowz
        avgdifferr = np.sqrt(differr_lowz**2. + differr_highz**2.)
        
        ax.errorbar(frvar.zHD,frvar.resid-frbase_matched.resid,yerr=np.sqrt(frvar.DLMAGERR**2.+frbase_matched.DLMAGERR**2.),fmt='.',color='k')
        ax.axhline(0.0,color='0.6',lw=2)
        ax.set_ylabel('$\Delta\mu$ (mag)')
        ax.text(
            0.5,0.73,f"{label}\n$\Delta\mu(z > 0.1) - \Delta\mu(z < 0.1) = {-1*avgdiff:.3f}\pm{avgdifferr:.3f}$",
            ha='center',va='bottom',
                transform=ax.transAxes,bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.8,'pad':0})
        ax.set_ylim([-0.7,0.7])

        axhist.xaxis.set_ticks([])
        axhist.set_ylabel('Hubble Resid.\n(mag)',labelpad=0)
        axhist.set_ylim([-0.7,0.7])

        mubins = np.linspace(-0.7,0.7,14)
        axhist.hist(frbase_matched.resid,bins=mubins,color='k',histtype='stepfilled',orientation='horizontal')
        axhist.hist(frvar.resid,bins=mubins,color='C0',histtype='stepfilled',orientation='horizontal')
        axhist.text(0.5,0.9,f"RMS={np.std(frbase_matched.resid):.3f}",
                    transform=axhist.transAxes,color='k',ha='center')
        axhist.text(0.5,0.8,f"RMS={np.std(frvar.resid):.3f}",
                    transform=axhist.transAxes,color='C0',ha='center')
        #import pdb; pdb.set_trace()
    ax4.set_xlabel('$z_{CMB}$',fontsize=15)

    # median low-z redshift: 0.02328
    # median high-z redshift: 0.42968
    # bins zero and 15 to look at median biascor difference
    # difference in median NIR biascor: -0.10485837220916425
    # difference in median opt+NIR biascor: -0.08549901874828461
    # difference in median Pantheon (SALT2) biascor: -0.03224136396116556

    #(Pdb) np.median(frvar.DLMAG[frvar.zHD > 0.15]-cosmo.mu(frvar.zHD[frvar.zHD > 0.15])) - np.median(frvar.DLMAG[frvar.zHD < 0.15]-cosmo.mu(frvar.zHD[frvar.zHD < 0.15]))
    #0.04364241308697814
    #(Pdb) np.median(frvar.DLMAG[frvar.zHD > 0.15]-cosmo.mu(frvar.zHD[frvar.zHD > 0.15])) - np.median(frvar.DLMAG[frvar.zHD < 0.15]-cosmo.mu(frvar.zHD[frvar.zHD < 0.15]))
    #0.09575035557295308
    
    import pdb; pdb.set_trace()

def table():

    frbase = txtobj(_nirfile,fitresheader=True)
    frbase.resid = frbase.DLMAG - cosmo.mu(frbase.zHD)
    mass_step_approx = np.average(frbase.resid[frbase.HOST_LOGMASS > 10],weights=1/frbase.DLMAGERR[frbase.HOST_LOGMASS > 10]**2.)-\
        np.average(frbase.resid[frbase.HOST_LOGMASS < 10],weights=1/frbase.DLMAGERR[frbase.HOST_LOGMASS < 10]**2.)
    print(mass_step_approx)
    frbase.resid[frbase.HOST_LOGMASS > 10] += np.abs(mass_step_approx)/2.
    frbase.resid[frbase.HOST_LOGMASS < 10] -= np.abs(mass_step_approx)/2.

    
    frbase.resid -= np.median(frbase.resid)
    frs = getmu.mkcuts(copy.deepcopy(_frs))
    #frs = copy.deepcopy(_frs)
    
    for frvar,label in zip(
            [_fropt,_frrv,_frs,_frnirmod],
            ['Optical+NIR ($R_V = 3.1$)',
             'Optical+NIR ($R_V = 2.0$)',
             'Optical with SALT2',
             'Optical+NIR $s_{BV}$, NIR dist.']):
        
        if 'SALT2' in label:
            frvar = getmu.getmu(frvar)
            frvar = getmu.mkcuts(frvar)
            frvar.resid = frvar.mures
            frvar.DLMAGERR = frvar.muerr
        
        frbase_matched = copy.deepcopy(frbase)
        frbase_matched.match_to_col('CID',frs.CID)
        frbase_matched.match_to_col('CID',frvar.CID)
        frvar.match_to_col('CID',frbase_matched.CID)

        if 'resid' not in frvar.__dict__.keys():
            frvar.resid = frvar.DLMAG - cosmo.mu(frvar.zHD)

        mass_step_approx = np.average(frvar.resid[frvar.HOST_LOGMASS > 10],weights=1/frvar.DLMAGERR[frvar.HOST_LOGMASS > 10]**2.)-\
            np.average(frvar.resid[frvar.HOST_LOGMASS < 10],weights=1/frvar.DLMAGERR[frvar.HOST_LOGMASS < 10]**2.)
        print(mass_step_approx)
            
        frvar.resid[frvar.HOST_LOGMASS > 10] += np.abs(mass_step_approx)/2.
        frvar.resid[frvar.HOST_LOGMASS < 10] -= np.abs(mass_step_approx)/2.
            
        frvar.resid -= np.median(frvar.resid)
        #if 'SALT2' in label: frvar.resid[frvar.zHD > 0.1] -= 0.073
        #else: frvar.resid[frvar.zHD > 0.1] -= 0.019

        diff_lowz,differr_lowz = weighted_avg_and_err(
            frvar.resid[frvar.zHD < 0.1]-frbase_matched.resid[frvar.zHD < 0.1],
            1/(frvar.DLMAGERR[frvar.zHD < 0.1]**2.+frbase_matched.DLMAGERR[frvar.zHD < 0.1]**2.))
        diff_highz,differr_highz = weighted_avg_and_err(
            frvar.resid[frvar.zHD > 0.1]-frbase_matched.resid[frvar.zHD > 0.1],
            1/(frvar.DLMAGERR[frvar.zHD > 0.1]**2.+frbase_matched.DLMAGERR[frvar.zHD > 0.1]**2.))
        avgdiff = diff_highz - diff_lowz
        avgdifferr = np.sqrt(differr_lowz**2. + differr_highz**2.)

        print(f"{label}&{-1*avgdiff:.3f}\pm{avgdifferr:.3f}$\\\\")
        
       # ax.text(
       #     0.5,0.73,f"{label}\n$\Delta\mu(z > 0.1) - \Delta\mu(z < 0.1) = {-1*avgdiff:.3f}\pm{avgdifferr:.3f}$",
       #     ha='center',va='bottom',
       #         transform=ax.transAxes,bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.8,'pad':0})


    # median low-z redshift: 0.02328
    # median high-z redshift: 0.42968
    # bins zero and 15 to look at median biascor difference
    # difference in median NIR biascor: -0.10485837220916425
    # difference in median opt+NIR biascor: -0.08549901874828461
    # difference in median Pantheon (SALT2) biascor: -0.03224136396116556

    #(Pdb) np.median(frvar.DLMAG[frvar.zHD > 0.15]-cosmo.mu(frvar.zHD[frvar.zHD > 0.15])) - np.median(frvar.DLMAG[frvar.zHD < 0.15]-cosmo.mu(frvar.zHD[frvar.zHD < 0.15]))
    #0.04364241308697814
    #(Pdb) np.median(frvar.DLMAG[frvar.zHD > 0.15]-cosmo.mu(frvar.zHD[frvar.zHD > 0.15])) - np.median(frvar.DLMAG[frvar.zHD < 0.15]-cosmo.mu(frvar.zHD[frvar.zHD < 0.15]))
    #0.09575035557295308
    
    import pdb; pdb.set_trace()

    
def newfig_hist():
    plt.rcParams['figure.figsize'] = (14,3)

    fig = plt.figure()

    
    plt.subplots_adjust(wspace=0,bottom=0.2,left=0.02,right=0.98)
    ax1hist = plt.subplot(151)
    ax2hist = plt.subplot(152)
    ax3hist = plt.subplot(153)
    ax4hist = plt.subplot(154)
    ax5hist = plt.subplot(155)
    #ax6hist = plt.subplot(166)
    
    frbase = txtobj(_nirfile,fitresheader=True)
    frbase.resid = frbase.DLMAG - cosmo.mu(frbase.zHD)
    mass_step_approx = np.average(frbase.resid[frbase.HOST_LOGMASS > 10],weights=1/frbase.DLMAGERR[frbase.HOST_LOGMASS > 10]**2.)-\
        np.average(frbase.resid[frbase.HOST_LOGMASS < 10],weights=1/frbase.DLMAGERR[frbase.HOST_LOGMASS < 10]**2.)
    print(mass_step_approx)
    frbase.resid[frbase.HOST_LOGMASS > 10] += np.abs(mass_step_approx)/2.
    frbase.resid[frbase.HOST_LOGMASS < 10] -= np.abs(mass_step_approx)/2.

    
    frbase.resid -= np.median(frbase.resid)
    #frs = getmu.mkcuts(copy.deepcopy(_frs))
    
    for axhist,frvar,label in zip(
            [ax1hist,ax2hist,ax3hist,ax4hist,ax5hist],
            [frbase,_fropt,_frrv,_frs,_frnirmod],            
            ['Baseline (NIR-only)',
             'Optical+NIR ($R_V = 1.5$)',
             'Optical+NIR ($R_V = 3.1$)',
             'Optical with SALT3',
             #'Optical with SALT3+cuts',
             'Optical+NIR $s_{BV}$, NIR dist.']):

        axhist.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        
        if 'SALT3' in label:
            frvar = getmu.getmu(frvar)
            if 'SALT3' in label:
                frvar = getmu.mkcuts(frvar,fitprobmin=0,x1errmin=2.0)

                frt = txtobj('output/fit_optical/CSP_RAISIN_SALT2_fixmjd.FITRES.TEXT',fitresheader=True)
                frt = getmu.getmu(frt)
                frt.resid = frt.mures
                frt.DLMAGERR = frt.muerr
                #import pdb; pdb.set_trace()
                for k in frvar.__dict__.keys():
                    frvar.__dict__[k] = np.append(frvar.__dict__[k],frt.__dict__[k][frt.CID == '2006is'])
                for k in frvar.__dict__.keys():
                    frvar.__dict__[k] = np.append(frvar.__dict__[k],frt.__dict__[k][frt.CID == '2009al'])

                #import pdb; pdb.set_trace()

            print(len(frvar.CID))
            frvar.resid = frvar.mures
            frvar.DLMAGERR = frvar.muerr
            #for i in frbase.CID:
            #    if i not in frvar.CID:
            #        print(i)
            #import pdb; pdb.set_trace()            
        #frbase_matched = copy.deepcopy(frbase)
        #frbase_matched.match_to_col('CID',frs.CID)
        #frbase_matched.match_to_col('CID',frvar.CID)
        #frvar.match_to_col('CID',frbase_matched.CID)

        if 'resid' not in frvar.__dict__.keys():
            frvar.resid = frvar.DLMAG - cosmo.mu(frvar.zHD)

        mass_step_approx = np.average(frvar.resid[frvar.HOST_LOGMASS > 10],weights=1/frvar.DLMAGERR[frvar.HOST_LOGMASS > 10]**2.)-\
            np.average(frvar.resid[frvar.HOST_LOGMASS < 10],weights=1/frvar.DLMAGERR[frvar.HOST_LOGMASS < 10]**2.)
        print(mass_step_approx)
            
        frvar.resid[frvar.HOST_LOGMASS > 10] += np.abs(mass_step_approx)/2.
        frvar.resid[frvar.HOST_LOGMASS < 10] -= np.abs(mass_step_approx)/2.

        frvar.resid -= np.median(frvar.resid)

        #diff_lowz,differr_lowz = weighted_avg_and_err(
        #    frvar.resid[frvar.zHD < 0.1]-frbase_matched.resid[frvar.zHD < 0.1],
        #    1/(frvar.DLMAGERR[frvar.zHD < 0.1]**2.+frbase_matched.DLMAGERR[frvar.zHD < 0.1]**2.))
        #diff_highz,differr_highz = weighted_avg_and_err(
        #    frvar.resid[frvar.zHD > 0.1]-frbase_matched.resid[frvar.zHD > 0.1],
        #    1/(frvar.DLMAGERR[frvar.zHD > 0.1]**2.+frbase_matched.DLMAGERR[frvar.zHD > 0.1]**2.))
        #avgdiff = diff_highz - diff_lowz
        #avgdifferr = np.sqrt(differr_lowz**2. + differr_highz**2.)
        
        axhist.yaxis.set_ticks([])

        axhist.set_xlim([-0.7,0.7])

        mubins = np.linspace(-0.7,0.7,14)
        #axhist.hist(frbase_matched.resid,bins=mubins,color='k',histtype='stepfilled')
        axhist.hist(frvar.resid,bins=mubins,color='0.5',histtype='bar',ec='0.0')
        axhist.hist(frvar.resid[frvar.zHD < 0.1],bins=mubins,color='b',alpha=0.5,hatch='\\',ec='0.0')
        axhist.hist(frvar.resid[frvar.zHD > 0.1],bins=mubins,color='r',alpha=0.5,hatch='//',ec='0.0')
        axhist.text(0.03,0.9,f"full sample RMS={np.std(frvar.resid):.3f}",fontsize=12,
                    transform=axhist.transAxes,color='0.2',ha='left',
                    bbox={'alpha':0.5,'facecolor':'1.0','edgecolor':'1.0','pad':0})
        axhist.text(0.03,0.8,fr"low-$z$ RMS={np.std(frvar.resid[frvar.zHD < 0.1]):.3f}",fontsize=12,
                    transform=axhist.transAxes,color='b',ha='left',
                    bbox={'alpha':0.5,'facecolor':'1.0','edgecolor':'1.0','pad':0})
        axhist.text(0.03,0.7,fr"RAISIN RMS={np.std(frvar.resid[frvar.zHD > 0.1]):.3f}",fontsize=12,
                    transform=axhist.transAxes,color='r',ha='left',
                    bbox={'alpha':0.5,'facecolor':'1.0','edgecolor':'1.0','pad':1})


        axhist.set_title(label)
        axhist.set_ylim([0,40])

    ax3hist.set_xlabel('Hubble Residual (mag)',fontsize=15)
    #xxl = ax4hist.set_xlabel('Hubble Resid. (mag)',fontsize=15)
    #xxl.set_position((xxl.get_position()[1],0.5))
    #xxl.set_horizontalalignment('center')    

    import pdb; pdb.set_trace()

    
def weighted_avg_and_err(values, weights):
    """
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, np.sqrt(variance/len(values)))

if __name__ == "__main__":
    #main()
    #newfig()
    newfig_hist()
    #table()
