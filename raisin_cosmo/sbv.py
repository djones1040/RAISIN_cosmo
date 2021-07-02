#!/usr/bin/env python
# D. Jones - 6/29/21

import pylab as plt
plt.rcParams['figure.figsize'] = (11,4)
plt.ion()
import numpy as np
from txtobj import txtobj
from scipy.optimize import minimize
from scipy.stats import binned_statistic
from scipy.odr import *

def neglnlikefunc(x,mu_i=None,sigma_i=None,sigma=None,survey=None,z=None):

    iCSP = survey == 'CSP'
    iPS1_midz = (survey == 'PS1') & (z < 0.4306)
    iPS1_highz = (survey == 'PS1') & (z >= 0.4306)
    iDES_midz = (survey == 'DES') & (z < 0.4306)
    iDES_highz = (survey == 'DES') & (z >= 0.4306)
    
    mu_lowz,mu_midz,mu_highz = x[0],x[1],x[2]
    #sigint_csp,sigint_ps1,sigint_des = x[3],x[4],x[5]
    sigint_csp,sigint_ps1,sigint_des = 0.15,0.15,0.15
    
    # sigint split by sample, but distance split by redshift
    # each one with a low-mass and high-mass component
    loglike_csp = -np.sum(-(mu_i[iCSP]-mu_lowz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                          np.log(1/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.))))
            
    loglike_ps1_midz = -np.sum(-(mu_i[iPS1_midz]-mu_midz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                               np.log(1//(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.))))
    
    loglike_ps1_highz = -np.sum(-(mu_i[iPS1_highz]-mu_highz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                np.log(1/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.))))
            
    loglike_des_midz = -np.sum(-(mu_i[iDES_midz]-mu_midz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                               np.log(1/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.))))
            
    loglike_des_highz = -np.sum(-(mu_i[iDES_highz]-mu_highz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                np.log(1/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.))))

    return loglike_csp + loglike_ps1_midz + loglike_ps1_highz + loglike_des_midz + loglike_des_highz

def errfnc(x):
    return(np.std(x)/np.sqrt(len(x)))

def main():

    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)

    fr = txtobj('st_av_corr_mags.txt')
    fr.SURVEY = np.array(['-99']*len(fr.CID))
    for j,i in enumerate(fr.CID):
        if fr.zHD[j] < 0.1: fr.SURVEY[j] = 'CSP'
        elif fr.CID[j].startswith('PS'): fr.SURVEY[j] = 'PS1'
        else: fr.SURVEY[j] = 'DES'

    # HR vs sBV (no sBV correction)
    # need to fit for and offset

    md = minimize(neglnlikefunc,(0,0.01,0.02), #0.09,0.1,0.11),
                  args=(fr.resid2,fr.DLMAGERR,None,fr.SURVEY,fr.zHD))
    fr.resid2[fr.SURVEY == 'CSP'] -= md.x[0]
    fr.resid2[(fr.zHD > 0.1) & (fr.zHD < 0.4306)] -= md.x[1]
    fr.resid2[(fr.zHD > 0.4306)] -= md.x[2]

    #ax1.errorbar()
    ax1.errorbar(fr.STRETCH,fr.resid2,xerr=fr.STRETCHERR,yerr=fr.DLMAGERR2,fmt='o',color='0.4')
    ax1.axhline(0,color='k')

    sbins = np.linspace(0.8,1.3,10)
    residbinned = binned_statistic(fr.STRETCH,fr.resid2,bins=sbins,statistic='median').statistic
    residbinnederr = binned_statistic(fr.STRETCH,fr.resid2,bins=sbins,statistic=errfnc).statistic
    ax1.errorbar((sbins[1:]+sbins[:-1])/2.,residbinned,yerr=residbinnederr,fmt='o-',color='C1')

    A = np.vander(fr.STRETCH, 2)
    C = np.diag(fr.DLMAGERR2**2.+0.15**2.)
    ATA = np.dot(A.T, A / (fr.DLMAGERR2**2.+0.15**2.)[:, None])
    cov = np.linalg.inv(ATA)
    w = np.linalg.solve(ATA, np.dot(A.T, fr.resid2 / (fr.DLMAGERR2**2.+0.15**2.)))
    print("Least-squares estimates:")
    print("m = {0:.3f} ± {1:.3f}".format(w[0], np.sqrt(cov[0, 0])))
    print("b = {0:.3f} ± {1:.3f}".format(w[1], np.sqrt(cov[1, 1])))

    def objective(x, a, b):
        return a * x + b
    s = np.linspace(0.8,1.3,100)
    resid_line = objective(s, w[0], w[1])
    resid_upper = objective(s, w[0]+np.sqrt(cov[0, 0]), w[1])
    resid_lower = objective(s, w[0]-np.sqrt(cov[0, 0]), w[1])
    ax1.plot(s,resid_line,color='r',ls='--',label=f'slope = {w[0]:.3f} ± {np.sqrt(cov[0, 0]):.3f}')
    #ax1.fill_between(s,resid_lower,resid_upper,color='0.7',alpha=0.5)
    ax1.legend()

    def quad_func(p, x):
        m, c = p
        return m*x + c
    quad_model = Model(quad_func)
    data = RealData(fr.STRETCH, fr.resid, sx=fr.STRETCHERR, sy=np.sqrt(fr.DLMAGERR**2.+0.15**2.))
    odr = ODR(data, quad_model, beta0=[0., 1.])
    out = odr.run()
    
    md = minimize(neglnlikefunc,(0,0.01,0.02), #0.09,0.1,0.11),
                  args=(fr.resid,fr.DLMAGERR,None,fr.SURVEY,fr.zHD))
    fr.resid[fr.SURVEY == 'CSP'] -= md.x[0]
    fr.resid[(fr.zHD > 0.1) & (fr.zHD < 0.4306)] -= md.x[1]
    fr.resid[(fr.zHD > 0.4306)] -= md.x[2]
    #fr.resid -= np.average(fr.resid,weights=1/fr.DLMAGERR**2.)
        
    # HR vs sBV (yes sBV correction)
    ax2.errorbar(fr.STRETCH,fr.resid,xerr=fr.STRETCHERR,yerr=fr.DLMAGERR,fmt='o',color='0.4')
    ax2.axhline(0,color='k')

    sbins = np.linspace(0.8,1.3,10)
    residbinned = binned_statistic(fr.STRETCH,fr.resid,bins=sbins,statistic='median').statistic
    residbinnederr = binned_statistic(fr.STRETCH,fr.resid,bins=sbins,statistic=errfnc).statistic
    ax2.errorbar((sbins[1:]+sbins[:-1])/2.,residbinned,yerr=residbinnederr,fmt='o-',color='C1')


    A = np.vander(fr.STRETCH, 2)
    C = np.diag(fr.DLMAGERR**2.+0.15**2.)
    ATA = np.dot(A.T, A / (fr.DLMAGERR ** 2+0.15**2.)[:, None])
    cov = np.linalg.inv(ATA)
    w = np.linalg.solve(ATA, np.dot(A.T, fr.resid / (fr.DLMAGERR ** 2+0.15**2.)))
    print("Least-squares estimates:")
    print("m = {0:.3f} ± {1:.3f}".format(w[0], np.sqrt(cov[0, 0])))
    print("b = {0:.3f} ± {1:.3f}".format(w[1], np.sqrt(cov[1, 1])))

    def objective(x, a, b):
        return a * x + b
    s = np.linspace(0.8,1.3,100)
    resid_line = objective(s, w[0], w[1])
    resid_upper = objective(s, w[0]+np.sqrt(cov[0, 0]), w[1])
    resid_lower = objective(s, w[0]-np.sqrt(cov[0, 0]), w[1])
    ax2.plot(s,resid_line,color='r',ls='--',label=f'slope = {w[0]:.3f} ± {np.sqrt(cov[0, 0]):.3f}')
    #ax2.fill_between(s,resid_lower,resid_upper,color='0.7',alpha=0.5)
    ax1.legend(prop={'size':12})
    ax2.legend(prop={'size':12})

    for ax in [ax1,ax2]:
        ax.set_xlabel('$s_{BV}$',fontsize=15,labelpad=0)
        ax.set_ylabel('Hubble Residual',fontsize=15)
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    ax1.set_title('Without $s_{BV}$ correction')
    ax2.set_title('With $s_{BV}$ correction')
    
    import pdb; pdb.set_trace()
    
if __name__ == "__main__":
    main()
