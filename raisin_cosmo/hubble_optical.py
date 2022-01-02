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
import scipy.stats
from scipy.optimize import minimize

_nirfile = 'output/cosmo_fitres/RAISIN_combined_FITOPT000_nobiascor.FITRES'
_csprvfile = 'output/fit_optical/CSP_RAISIN_optnir_MWRV.FITRES.TEXT'
_ps1rvfile = 'output/fit_optical/PS1_RAISIN_optnir_MWRV.FITRES.TEXT'
_desrvfile = 'output/fit_optical/DES_RAISIN_optnir_MWRV.FITRES.TEXT'
_cspoptfile = 'output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT'
_ps1optfile = 'output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT'
_desoptfile = 'output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT'
_cspoptonlyfile = 'output/fit_optical/CSP_RAISIN_optical.FITRES.TEXT'
_ps1optonlyfile = 'output/fit_optical/PS1_RAISIN_optical.FITRES.TEXT'
_desoptonlyfile = 'output/fit_optical/DES_RAISIN_optical.FITRES.TEXT'

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
froptonlycsp = txtobj(_cspoptonlyfile,fitresheader=True)
froptonlyps1 = txtobj(_ps1optonlyfile,fitresheader=True)
froptonlydes = txtobj(_desoptonlyfile,fitresheader=True)

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
    if k != 'RV' and k != 'RVERR':
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

# here's optical only
idxcsp,idxps1,idxdes = np.array([],dtype=int),np.array([],dtype=int),np.array([],dtype=int)
for j,i in enumerate(frnir.CID):
    if i not in _goodcids: continue
    if i in froptonlycsp.CID:
        idxcsp = np.append(idxcsp,np.where(froptonlycsp.CID == i)[0][0])
    if i in froptonlyps1.CID:
        idxps1 = np.append(idxps1,np.where(froptonlyps1.CID == i)[0][0])
    if i in froptonlydes.CID:
        idxdes = np.append(idxdes,np.where(froptonlydes.CID == i)[0][0])

if 'PKMJDINI' in froptonlycsp.__dict__.keys():
    del froptonlycsp.__dict__['PKMJDINI']
    del froptonlycsp.__dict__['COV_STRETCH_AV']
    del froptonlycsp.__dict__['COV_AV_DLMAG']
for k in froptonlycsp.__dict__.keys():
    froptonlycsp.__dict__[k] = np.concatenate(
        (froptonlycsp.__dict__[k][idxcsp],froptonlyps1.__dict__[k][idxps1],
        froptonlydes.__dict__[k][idxdes]))
_froptonly = froptonlycsp


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

def apply_all_cuts(fropt):

    # AV
    iGoodAV = np.zeros(len(fropt.CID),dtype=bool)
    for j,i in enumerate(fropt.CID):
        if i in fropt.CID and fropt.AV[fropt.CID == i] < 0.3*fropt.RV[fropt.CID == i]:
            iGoodAV[j] = True

    # reasonable stretch
    iGoodSt = np.zeros(len(fropt.CID),dtype=bool)
    for j,i in enumerate(fropt.CID):
        if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > 0.75 and fropt.STRETCH[fropt.CID == i] < 1.185 and \
           fropt.STRETCHERR[fropt.CID == i] < 0.3: #175:
            iGoodSt[j] = True

    #iGoodSt = (fropt.STRETCH > 0.8) & (fropt.STRETCH < 1.3)

    for k in fropt.__dict__.keys():
        fropt.__dict__[k] = fropt.__dict__[k][iGoodAV & iGoodSt]

    return fropt


def lnlikefunc(x,p_iae=None,mu_i=None,sigma_i=None,sigma=None,z=None,survey=None):

    if sigma or sigma == 0.0:
        # fix the dispersion
        x[2] = sigma; x[3] = sigma

    p_iae[np.where(p_iae == 0)] == 1e-4
    return -np.sum(np.logaddexp(-(mu_i-x[0])**2./(2.0*(sigma_i**2.+x[2]**2.)) +\
                np.log((1-0.01*p_iae)/(np.sqrt(2*np.pi)*np.sqrt(x[2]**2.+sigma_i**2.))),
                -(mu_i-x[1])**2./(2.0*(sigma_i**2.+x[3]**2.)) +\
                np.log((0.01*p_iae)/(np.sqrt(2*np.pi)*np.sqrt(x[3]**2.+sigma_i**2.)))))

def neglnlikefunc(x,p_lm=None,mu_i=None,sigma_i=None,sigma=None,survey=None,z=None):

    iCSP = survey == 5
    iPS1_midz = (survey == 15) & (z < 0.4306)
    iPS1_highz = (survey == 15) & (z >= 0.4306)
    iDES_midz = (survey == 10) & (z < 0.4306)
    iDES_highz = (survey == 10) & (z >= 0.4306)
    
    mu_lowz,mu_midz,mu_highz = x[0],x[1],x[2]
    #mu_midz = mu_highz
    #mu_lowz = mu_highz
    sigint_csp,sigint_ps1,sigint_des = x[3],x[4],x[5]
    mass_step = x[6]
    
    # sigint split by sample, but distance split by redshift
    # each one with a low-mass and high-mass component
    loglike_csp = -np.sum(np.logaddexp(-(mu_i[iCSP]+mass_step-mu_lowz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                                       np.log((1-0.01*p_lm[iCSP])/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.))),
                                       -(mu_i[iCSP]-mu_lowz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                                       np.log((0.01*p_lm[iCSP])/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.)))))
            
    loglike_ps1_midz = -np.sum(np.logaddexp(-(mu_i[iPS1_midz]+mass_step-mu_midz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                                            np.log((1-0.01*p_lm[iPS1_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.))),
                                            -(mu_i[iPS1_midz]-mu_midz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                                            np.log((0.01*p_lm[iPS1_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.)))))
    
    loglike_ps1_highz = -np.sum(np.logaddexp(-(mu_i[iPS1_highz]+mass_step-mu_highz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                             np.log((1-0.01*p_lm[iPS1_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.))),
                                             -(mu_i[iPS1_highz]-mu_highz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                             np.log((0.01*p_lm[iPS1_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.)))))
            
    loglike_des_midz = -np.sum(np.logaddexp(-(mu_i[iDES_midz]+mass_step-mu_midz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                                            np.log((1-0.01*p_lm[iDES_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.))),
                                            -(mu_i[iDES_midz]-mu_midz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                                            np.log((0.01*p_lm[iDES_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.)))))
            
    loglike_des_highz = -np.sum(np.logaddexp(-(mu_i[iDES_highz]+mass_step-mu_highz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                             np.log((1-0.01*p_lm[iDES_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.))),
                                             -(mu_i[iDES_highz]-mu_highz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                             np.log((0.01*p_lm[iDES_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.)))))

    return loglike_csp + loglike_ps1_midz + loglike_ps1_highz + loglike_des_midz + loglike_des_highz


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
    #mass_step_approx = np.average(frbase.resid[frbase.HOST_LOGMASS > 10],weights=1/frbase.DLMAGERR[frbase.HOST_LOGMASS > 10]**2.)-\
    #    np.average(frbase.resid[frbase.HOST_LOGMASS < 10],weights=1/(frbase.DLMAGERR[frbase.HOST_LOGMASS < 10]**2.+0.1**2.))
    #print(mass_step_approx)
    #frbase.resid[frbase.HOST_LOGMASS > 10] += np.abs(mass_step_approx)/2.
    #frbase.resid[frbase.HOST_LOGMASS < 10] -= np.abs(mass_step_approx)/2.

    #import pdb; pdb.set_trace()
    #frbase.resid -= np.median(frbase.resid)
    frs = getmu.mkcuts(copy.deepcopy(_frs),fitprobmin=0)
    #frs = copy.deepcopy(_frs)
    
    for frvar,label in zip(
            [_fropt,_froptonly,_frrv,_frs,_frnirmod],
            ['Optical+NIR ($R_V = 1.5$)',
             'Optical Only ($R_V = 1.5$)',
             'Optical+NIR ($R_V = 3.1$)',
             'Optical with SALT2',
             'Optical+NIR $s_{BV}$, NIR dist.']):
        
        if 'SALT2' in label:
            frvar = getmu.getmu(frvar)
            frvar = getmu.mkcuts(frvar)
            frvar.resid = frvar.mures
            frvar.DLMAGERR = frvar.muerr
        #if 'Optical+NIR ($R' in label:
        #    frvar = apply_all_cuts(frvar)
            
        frbase_matched = copy.deepcopy(frbase)
        if 'SALT2' in label: frbase_matched.match_to_col('CID',frs.CID)
        frbase_matched.match_to_col('CID',frvar.CID)
        frvar.match_to_col('CID',frbase_matched.CID)

        if 'resid' not in frvar.__dict__.keys():
            frvar.resid = frvar.DLMAG - cosmo.mu(frvar.zHD)

        frbase_matched.HOST_LOGMASS_ERR = np.sqrt(frbase_matched.HOST_LOGMASS_ERR**2. + 0.02**2.)
        frbase_matched.p_hm = np.zeros(len(frbase_matched.CID))
        boundary = 10
        for i in range(len(frbase_matched.CID)):
            frbase_matched.p_hm[i] = scipy.stats.norm.cdf(
                boundary,float(frbase_matched.HOST_LOGMASS[i]),
                float(frbase_matched.HOST_LOGMASS_ERR[i]))*100.
        md = minimize(neglnlikefunc,(0.0,0.0,0.1,0.10,0.01,0.02,0.09,0.1,0.11,0.1),
                      args=(frbase_matched.p_hm,frbase_matched.resid,
                            frbase_matched.DLMAGERR,None,frbase_matched.IDSURVEY,frbase_matched.zHD))
        frbase_matched.resid[frbase_matched.HOST_LOGMASS > 10] += md.x[6]/2.
        frbase_matched.resid[frbase_matched.HOST_LOGMASS < 10] -= md.x[6]/2.
#        import pdb; pdb.set_trace()
            
        frvar.HOST_LOGMASS_ERR = np.sqrt(frbase_matched.HOST_LOGMASS_ERR**2. + 0.02**2.)
        frvar.IDSURVEY = frbase_matched.IDSURVEY
        frvar.p_hm = np.zeros(len(frvar.CID))
        boundary = 10
        for i in range(len(frvar.CID)):
            frvar.p_hm[i] = scipy.stats.norm.cdf(
                boundary,float(frvar.HOST_LOGMASS[i]),
                float(frvar.HOST_LOGMASS_ERR[i]))*100.
        md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
                      args=(frvar.p_hm,frvar.resid,frvar.DLMAGERR,None,frvar.IDSURVEY,frvar.zHD))
        frvar.resid[frvar.HOST_LOGMASS > 10] += md.x[6]/2.
        frvar.resid[frvar.HOST_LOGMASS < 10] -= md.x[6]/2.


        diff_lowz,differr_lowz = weighted_avg_and_err(
            frvar.resid[frvar.zHD < 0.1]-frbase_matched.resid[frvar.zHD < 0.1],
            1/(frvar.DLMAGERR[frvar.zHD < 0.1]**2.+frbase_matched.DLMAGERR[frvar.zHD < 0.1]**2.+0.1**2.))

        diff_highz,differr_highz = weighted_avg_and_err(
            frvar.resid[frvar.zHD > 0.1]-frbase_matched.resid[frvar.zHD > 0.1],
            1/(frvar.DLMAGERR[frvar.zHD > 0.1]**2.+frbase_matched.DLMAGERR[frvar.zHD > 0.1]**2.+0.1**2.))

        avgdiff = diff_highz - diff_lowz
        avgdifferr = np.sqrt(differr_lowz**2. + differr_highz**2.)

        print(f"{label}&{-1*avgdiff:.3f}\pm{avgdifferr:.3f}$\\\\")

        #plt.errorbar(frbase_matched.zHD,frbase_matched.DLMAG-cosmo.mu(frbase_matched.zHD),
        #             yerr=frbase_matched.DLMAGERR,fmt='o')
        #plt.errorbar(frvar.zHD,frvar.DLMAG-cosmo.mu(frvar.zHD),yerr=frvar.DLMAGERR,fmt='o')
        #plt.axhline(0,color='k')
        import pdb; pdb.set_trace()
        #frvar_back = copy.deepcopy(frvar)
        plt.clf()
        #np.median(frvar.DLMAG[frvar.zHD < 0.1]-cosmo.mu(frvar.zHD[frvar.zHD < 0.1]))
        #np.median(frvar.DLMAG[frvar.zHD > 0.1]-cosmo.mu(frvar.zHD[frvar.zHD > 0.1]))
        #np.median(frbase_matched.DLMAG[frvar.zHD < 0.1]-cosmo.mu(frvar.zHD[frvar.zHD < 0.1]))
        #np.median(frbase_matched.DLMAG[frvar.zHD > 0.1]-cosmo.mu(frvar.zHD[frvar.zHD > 0.1]))
        
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

    
    plt.subplots_adjust(wspace=0,bottom=0.2,left=0.05,right=0.98)
    ax1hist = plt.subplot(151)
    ax2hist = plt.subplot(152)
    ax3hist = plt.subplot(153)
    ax4hist = plt.subplot(155)
    ax5hist = plt.subplot(154)
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
        #print(len(frvar.resid))
        mass_step_approx = np.average(frvar.resid[frvar.HOST_LOGMASS > 10],weights=1/frvar.DLMAGERR[frvar.HOST_LOGMASS > 10]**2.)-\
            np.average(frvar.resid[frvar.HOST_LOGMASS < 10],weights=1/frvar.DLMAGERR[frvar.HOST_LOGMASS < 10]**2.)
        print(mass_step_approx)
        print('hack')
        #frvar.resid[frvar.HOST_LOGMASS > 10] += np.abs(mass_step_approx)/2.
        #frvar.resid[frvar.HOST_LOGMASS < 10] -= np.abs(mass_step_approx)/2.

        #frvar.resid -= np.median(frvar.resid)

        #diff_lowz,differr_lowz = weighted_avg_and_err(
        #    frvar.resid[frvar.zHD < 0.1]-frbase_matched.resid[frvar.zHD < 0.1],
        #    1/(frvar.DLMAGERR[frvar.zHD < 0.1]**2.+frbase_matched.DLMAGERR[frvar.zHD < 0.1]**2.))
        #diff_highz,differr_highz = weighted_avg_and_err(
        #    frvar.resid[frvar.zHD > 0.1]-frbase_matched.resid[frvar.zHD > 0.1],
        #    1/(frvar.DLMAGERR[frvar.zHD > 0.1]**2.+frbase_matched.DLMAGERR[frvar.zHD > 0.1]**2.))
        #avgdiff = diff_highz - diff_lowz
        #avgdifferr = np.sqrt(differr_lowz**2. + differr_highz**2.)
        
        #axhist.yaxis.set_ticks([])
        if axhist != ax1hist: axhist.yaxis.set_ticklabels([])
        
        axhist.set_xlim([-0.7,0.7])

        mubins = np.linspace(-0.7,0.7,14)
        #axhist.hist(frbase_matched.resid,bins=mubins,color='k',histtype='stepfilled')
        axhist.hist(frvar.resid,bins=mubins,color='0.5',histtype='bar',ec='0.0')
        axhist.hist(frvar.resid[frvar.zHD < 0.1],bins=mubins,color='b',alpha=0.5,hatch='\\',ec='0.0')
        axhist.hist(frvar.resid[frvar.zHD > 0.1],bins=mubins,color='r',alpha=0.5,hatch='//',ec='0.0')
        axhist.text(0.03,0.9,f"full sample RMS",fontsize=12,
                    transform=axhist.transAxes,color='0.2',ha='left',
                    bbox={'alpha':0.5,'facecolor':'1.0','edgecolor':'1.0','pad':0})
        axhist.text(0.03,0.8,fr"low-$z$ RMS",fontsize=12,
                    transform=axhist.transAxes,color='b',ha='left',
                    bbox={'alpha':0.5,'facecolor':'1.0','edgecolor':'1.0','pad':0})
        axhist.text(0.03,0.7,fr"RAISIN RMS",fontsize=12,
                    transform=axhist.transAxes,color='r',ha='left',
                    bbox={'alpha':0.5,'facecolor':'1.0','edgecolor':'1.0','pad':1})

        axhist.text(0.97,0.9,f"{np.std(frvar.resid):.3f}",fontsize=12,
                    transform=axhist.transAxes,color='0.2',ha='right',
                    bbox={'alpha':0.5,'facecolor':'1.0','edgecolor':'1.0','pad':0})
        axhist.text(0.97,0.8,fr"{np.std(frvar.resid[frvar.zHD < 0.1]):.3f}",fontsize=12,
                    transform=axhist.transAxes,color='b',ha='right',
                    bbox={'alpha':0.5,'facecolor':'1.0','edgecolor':'1.0','pad':0})
        axhist.text(0.97,0.7,fr"{np.std(frvar.resid[frvar.zHD > 0.1]):.3f}",fontsize=12,
                    transform=axhist.transAxes,color='r',ha='right',
                    bbox={'alpha':0.5,'facecolor':'1.0','edgecolor':'1.0','pad':1})

        

        cidlist_test = ['PScA470110', 'PScA470240', 'PScB480464', 'PScB480794',
                        'PScC490037', 'PScD500100', 'PScD500301', 'PScF510457',
                        'PScF520062', 'PScF520188', 'PScG530251', 'PScH540087',
                        'PScH540118', 'PScJ440005', 'PScJ440236', 'PScJ550202',
                        'PScJ560027', 'PScJ560054', 'PScK450339', 'DES15C3odz',
                        'DES15E2mhy', 'DES15E2nlz', 'DES15E2uc', 'DES15X2kvt',
                        'DES15X2mey', 'DES15X2nkz', 'DES16C1cim', 'DES16C2cva',
                        'DES16C3cmy', 'DES16E2clk', 'DES16E2cqq', 'DES16S1agd',
                        'DES16S1bno', 'DES16S2afz', 'DES16X3cry', 'DES16X3zd',
                        'DES16E1dcx']
        dlmaglist_test = [41.31596, 41.50991, 40.03098, 41.30598, 41.61074, 40.86769,
       40.9813 , 42.59389, 41.30391, 40.6775 , 41.68742, 40.78568,
       42.11243, 40.87097, 41.77584, 41.59896, 41.44466, 41.9577 ,
       41.61716, 42.38261, 41.73083, 41.82508, 42.77946, 41.67113,
       42.8961 , 42.08552, 42.45433, 41.70474, 42.43624, 41.50161,
       41.92863, 42.04531, 42.18525, 41.96995, 42.76736, 42.02113,
       41.87013]
        zlist_test = [0.3456, 0.4306, 0.2212, 0.3342, 0.4231, 0.3106, 0.3254, 0.5027,
       0.308 , 0.2804, 0.4118, 0.275 , 0.4766, 0.3056, 0.4292, 0.4208,
       0.4392, 0.4804, 0.4096, 0.508 , 0.4391, 0.41  , 0.566 , 0.404 ,
       0.608 , 0.4688, 0.531 , 0.4029, 0.5564, 0.367 , 0.426 , 0.504 ,
       0.47  , 0.483 , 0.612 , 0.495 , 0.453 ]
        residlist_test = [ 0.00403203, -0.36402382, -0.17057959,  0.0789223 , -0.2179347 ,
       -0.17486724, -0.17840375,  0.31867153,  0.28246053, -0.10918929,
       -0.07164382,  0.04739309, -0.02409095, -0.13084958, -0.08970109,
       -0.21568655, -0.48028187, -0.19944762, -0.12814631,  0.08004832,
       -0.19352408,  0.07726725,  0.19393735, -0.03884689,  0.12218973,
       -0.00817877,  0.03611209,  0.00175591, -0.10438031,  0.03709405,
        0.08236778, -0.23664008,  0.08492032, -0.20122338, -0.02384703,
       -0.21388083, -0.13476052]

        #residlist_out = []
        #for i,c in enumerate(cidlist_test):
        #    if abs(frvar.DLMAG[frvar.zHD > 0.1][frvar.CID[frvar.zHD > 0.1] == c]-dlmaglist_test[i]) > 1e-4:
        #        print('hi1',c)
        #    if abs(frvar.zHD[frvar.zHD > 0.1][frvar.CID[frvar.zHD > 0.1] == c]-zlist_test[i]) > 1e-4:
        #        print('hi2',c)
        #    if abs(frvar.resid[frvar.zHD > 0.1][frvar.CID[frvar.zHD > 0.1] == c]-residlist_test[i]) > 1e-4:
        #        print('hi3',c)
        #    else:
        #        try: residlist_out += [frvar.resid[frvar.zHD > 0.1][frvar.CID[frvar.zHD > 0.1] == c][0]]
        #        except:
        #            import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()

        axhist.set_title(label)
        axhist.set_ylim([0,40])

    ax3hist.set_xlabel('Hubble Residual (mag)',fontsize=15)
    #xxl = ax4hist.set_xlabel('Hubble Resid. (mag)',fontsize=15)
    #xxl.set_position((xxl.get_position()[1],0.5))
    #xxl.set_horizontalalignment('center')    
    ax1hist.set_ylabel(r'N$_{\rm SNe}$',fontsize=15)
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
