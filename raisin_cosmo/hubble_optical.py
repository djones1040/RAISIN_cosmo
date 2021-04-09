#!/usr/bin/env python
# D. Jones - 1/27/21

from txtobj import txtobj
import numpy as np
import pylab as plt
plt.ion()
import cosmo
import scipy

_nirfile = 'output/cosmo_fitres/RAISIN_combined_FITOPT000_nobiascor.FITRES'
_cspoptfile = 'output/fit_optical/CSP_RAISIN_optnir_LOWRV.FITRES.TEXT'
_ps1optfile = 'output/fit_optical/PS1_RAISIN_optnir_LOWRV.FITRES.TEXT'
_desoptfile = 'output/fit_optical/DES_RAISIN_optnir_LOWRV.FITRES.TEXT'
_csprvfile = 'output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT'
_ps1rvfile = 'output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT'
_desrvfile = 'output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT'
#_csprvfile = 'output/fit_optical/CSP_RAISIN_SALT2.FITRES.TEXT'
#_ps1rvfile = 'output/fit_optical/PS1_RAISIN_SALT2.FITRES.TEXT'
#_desrvfile = 'output/fit_optical/DES_RAISIN_SALT2.FITRES.TEXT'

# NIR vs. optical distances
# NIR vs. optical with RV floated

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

    fropt.DLMAG -= np.median(fropt.DLMAG-cosmo.mu(fropt.zHD))
    frrv.DLMAG -= np.median(frrv.DLMAG-cosmo.mu(frrv.zHD))
    frnir.DLMAG -= np.median(frnir.DLMAG-cosmo.mu(frnir.zHD))    
    ax1.errorbar(fropt.zHD,fropt.DLMAG-cosmo.mu(fropt.zHD),yerr=fropt.DLMAGERR,fmt='*',color='b',label=f'Optical+NIR, $R_V = 2.0$, med. $\Delta\mu = {delmuopt:.3f}$')
    ax1.errorbar(frrv.zHD,frrv.DLMAG-cosmo.mu(frrv.zHD),yerr=frrv.DLMAGERR,fmt='D',color='lightblue',label=f'Optical+NIR, $R_V = 3.1$, med. $\Delta\mu = {delmurv:.3f}$')
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

    
    
if __name__ == "__main__":
    main()
