#!/usr/bin/env python

from txtobj import txtobj
import numpy as np
import cosmo
import matplotlib.pyplot as plt
import numpy.linalg as la
import os
plt.ion()
filtlist = 'BgriYJH'
#lowzpars = txtobj('../GPfit/lcparams_avelino19.txt')
#badcid = ['sn2006N','sn2005lu','sn2010iw','sn2004S','sn2008gb',
#          'sn2008gg','sn2007st','sn2009kk','sn2007as','sn2008C','sn2008fl',
#          'sn2005na','sn2007ai','sn2009ag','sn2008fr','sn2008af','sn2001cn',
#          'sn2000ca','sn2003hv','snf20080514-002','sn2000bh',
#          'sn2006ej','sn2007sr','sn2008fw']
# aaaaaand things with Bs,Vs,Rs,Is, see if removing those helps
#badcid += ['sn2003du','sn2010iw','sn2004S','sn2009kw', #maybe?
#          'sn2009an','sn2006cp','sn2007co','sn2011by',
#          'sn2011df','sn2008hm','sn2002dj','sn2001el',
#          'sn2006lf','sn2006D','sn1999ek','sn2009al',
#          'sn2009kk','sn1999ee','sn1998bu','sn2010dw',
#          'sn2006ac','sn2007qe','sn2008hs','sn2001ba',
#          'sn2001bt','sn2010ai','sn2008af','sn2001cn',
#          'sn2000ca','sn2000E','sn2003hv','sn2009bv',
#          'sn2001cz','sn2000bh','sn2011B','sn2011ao','sn2007cq','sn2010kg']

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, np.sqrt(variance))

def isPD(B):
    """Returns true when input is positive-definite, via Cholesky"""
    try:
        _ = la.cholesky(B)
        return True
    except la.LinAlgError:
        return False

def nearestPD(A):
    """Find the nearest positive-definite matrix to input

    A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
    credits [2].

    [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

    [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
    """

    B = (A + A.T) / 2
    _, s, V = la.svd(B)

    H = np.dot(V.T, np.dot(np.diag(s), V))

    A2 = (B + H) / 2

    A3 = (A2 + A2.T) / 2
    A3 = np.around(A3,decimals=3)

    if isPD(A3):
        return A3

    spacing = np.spacing(la.norm(A))
    # The above is different from [1]. It appears that MATLAB's `chol` Cholesky
    # decomposition will accept matrixes with exactly 0-eigenvalue, whereas
    # Numpy's will not. So where [1] uses `eps(mineig)` (where `eps` is Matlab
    # for `np.spacing`), we use the above definition. CAVEAT: our `spacing`
    # will be much larger than [1]'s `eps(mineig)`, since `mineig` is usually on
    # the order of 1e-16, and `eps(1e-16)` is on the order of 1e-34, whereas
    # `spacing` will, for Gaussian random matrixes of small dimension, be on
    # othe order of 1e-16. In practice, both ways converge, as the unit test
    # below suggests.
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(la.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1
        A3 = np.around(A3,decimals=3)

    return A3

def plothubble():

    fr = txtobj('opticalnir_maxmodel_fitresults_lowz.txt')
    plt.plot(fr.zcmb,fr.Jmu)
    
    import pdb; pdb.set_trace()
    
def main():
    
    plt.rcParams['figure.figsize'] = (7,7)
    fr = txtobj('opticalnir_maxmodel_fitresults_lowz.txt')
    fr.zcmb = np.zeros(len(fr.z))
    fr.ebv = np.zeros(len(fr.z))
    iGood = np.array([],dtype=int)
    for i in range(len(fr.SNID)):
        fr.zcmb[i] = lowzpars.z_cmb[lowzpars.snid == 'SN'+fr.SNID[i][2:]][0]
        fr.ebv[i] = lowzpars.EBVMW[lowzpars.snid == 'SN'+fr.SNID[i][2:]][0]
        if fr.SNID[i] not in badcid and fr.zcmb[i] > 0.015 and fr.st[i] < 1.3 and fr.st[i] > 0.7 and fr.EBV[i] < 0.3:
            iGood = np.append(iGood,i)
    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iGood]
            
    # need to shape and EBV-correct these
    for filt in filtlist:
        iGood = fr.__dict__['%smu'%filt] != -99
        zerr = 0.00083*5.0/np.log(10)*(1.0+fr.zcmb[iGood])/(fr.zcmb[iGood]*(1.0+fr.zcmb[iGood]/2.0))

        fr.__dict__['%smures'%filt] = np.array([-99.]*len(fr.__dict__['%smu'%filt]))
        fr.__dict__['%smures'%filt][iGood] = fr.__dict__['%smu'%filt][iGood] - cosmo.mu(fr.zcmb[iGood])
        fr.__dict__['%smures'%filt][iGood] -= np.average(fr.__dict__['%smures'%filt][iGood],weights=1/(fr.__dict__['%smuerr'%filt][iGood]**2.+zerr**2.))

        
    covmat = np.zeros([7,7])
    for i,flt1 in enumerate(filtlist):
        for j,flt2 in enumerate(filtlist):
            imures = (fr.__dict__['%smures'%flt1] != -99) & (fr.__dict__['%smures'%flt2] != -99)
            if len(np.where(imures)[0]) > 2:
                covmat[j,i] = np.sum(fr.__dict__['%smures'%flt1][imures]*fr.__dict__['%smures'%flt2][imures])/len(fr.__dict__['%smures'%flt1][imures])
                #import pdb; pdb.set_trace()
            else:
                covmat[j,i] = 0.0
            #if flt1 == 'Y' and flt2 == 'Y': import pdb; pdb.set_trace()


    print(np.around(np.sqrt(np.diag(covmat)),decimals=3))
    print(np.around(cov2corr(covmat),decimals=3))
    print(nearestPD(cov2corr(covmat)))
    
    ax = plt.axes([0.1,0.07,0.8,0.7])
    cax = plt.axes([0.1,0.82,0.8,0.1])
    cb = ax.imshow(covmat)
    clb = plt.colorbar(cb,cax=cax,orientation='horizontal')
    cax.set_xlabel('$\sigma_{\mu}^2$',labelpad=-100,fontsize=20)
    
    ax.xaxis.set_ticks([0,1,2,3,4,5,6])
    ax.yaxis.set_ticks([0,1,2,3,4,5,6])
    ax.xaxis.set_ticklabels(['B','g','r','i','Y','J','H'])
    ax.yaxis.set_ticklabels(['B','g','r','i','Y','J','H'])
    import pdb; pdb.set_trace()

_nmltmpl = """
  &SNLCINP
     PRIVATE_DATA_PATH = '$RAISIN_ROOT/cosmo/data/Photometry'

     VERSION_PHOTOMETRY = '<data_version>'
     KCOR_FILE         = '$RAISIN_ROOT/cosmo/kcor/<kcor>'

     NFIT_ITERATION = 3
     INTERP_OPT     = 1

     SNTABLE_LIST = 'FITRES LCPLOT(text:key)'
     TEXTFILE_PREFIX  = '<outfile>'
     
     LDMP_SNFAIL = T
     USE_MWCOR = F
     USE_MINOS = F

     H0_REF   = 70.0
     OLAM_REF =  0.70
     OMAT_REF =  0.30
     W0_REF   = -1.00

     SNCID_LIST    =  0
     CUTWIN_CID    =  0, 20000000
     SNCCID_LIST   =  <cidlist>
     SNCCID_IGNORE =  

     cutwin_redshift   = 0.001, 2.0
     cutwin_Nepoch    =  1

     RV_MWCOLORLAW = 3.1
     !OPT_MWCOLORLAW = 99
     OPT_MWEBV = 3
     MWEBV_SCALE = 1.00
     MWEBV_SHIFT = 0.0
     FUDGE_MAG_ERROR = 


     MAGOBS_SHIFT_PRIMARY = ' '
     EPCUT_SNRMIN = ''
     ABORT_ON_NOEPOCHS = F
     HEADER_OVERRIDE_FILE= '$RAISIN_ROOT/cosmo/vpec_baseline_raisin.list'

  &END
  &FITINP

     FITMODEL_NAME  = '$RAISIN_ROOT/cosmo/snoopy.B18'
     OPT_PRIOR_AV = 0
    
     PRIOR_MJDSIG        = 5.0
     PRIOR_LUMIPAR_RANGE = -5.0, 5.0
     INIVAL_GRIDSEARCH_COLOR = -1.0, 1.0, 0.05
     !INIVAL_SHAPE = <inival_st>
     !INISTP_SHAPE = 0.0
     !INIVAL_AV = <inival_av>
     !INISTP_AV = 0.0
     !INIVAL_PEAKMJD = <inival_pkmjd>
     !INISTP_PEAKMJD = 0.0
     !INIVAL_DLMAG = <inival_dlmag>
     !INISTP_DLMAG = 0.0

     INIVAL_RV = 1.518

     !OPT_COVAR = 1
     OPT_XTMW_ERR = 1
     OPT_COVAR_FLUX = 0
     TREST_REJECT  = -15.0, 45.0
     NGRID_PDF     = 0

     FUDGEALL_ITER1_MAXFRAC = 0.02
     FILTLIST_FIT = '<filtlist>'

  &END
"""

    
class dispmodel:

    def __init__(self):
        self.goodcids = np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)

    def run_snana(self,cidlist,version,kcor,filtlist,outfile,
                  inival_st=None,inival_av=None,inival_dlmag=None,
                  inival_pkmjd=None,clobber=False):

        if not clobber and os.path.exists(f"{outfile}.FITRES.TEXT"):
            return outfile
        
        nmltext = _nmltmpl.replace('<data_version>',version).\
            replace('<kcor>',kcor).\
            replace('<filtlist>',filtlist).\
            replace('<outfile>',outfile).\
            replace('<cidlist>',"'"+"','".join(cidlist)+"'")

        if inival_st is not None:
            nmltext = nmltext.replace('<inival_st>',f"{inival_st:.3f}").\
                replace('!INISTP_SHAPE','INISTP_SHAPE').\
                replace('!INIVAL_SHAPE','INIVAL_SHAPE').\
                replace('!INISTP_PEAKMJD','INISTP_PEAKMJD')
        else:
            nmltext = nmltext.replace('<inival_st>','1.0')

        if inival_av is not None:
            nmltext = nmltext.replace('<inival_av>',f"{inival_av:.3f}").\
                replace('!INISTP_AV','INISTP_AV').\
                replace('!INIVAL_AV','INIVAL_AV').\
                replace('!INISTP_PEAKMJD','INISTP_PEAKMJD')
        else:
            nmltext = nmltext.replace('<inival_av>','1.0')

        if inival_pkmjd is not None:
            nmltext = nmltext.replace('<inival_pkmjd>',f"{inival_pkmjd:.3f}").\
                replace('!INIVAL_PEAKMJD','INIVAL_PEAKMJD')
        else:
            pass

            
        if inival_dlmag is not None:
            nmltext = nmltext.replace('<inival_dlmag>',f"{inival_dlmag:.3f}").\
                replace('!INISTP_DLMAG','INISTP_DLMAG').\
                replace('!INIVAL_DLMAG','INIVAL_DLMAG').\
                replace('!INISTP_PEAKMJD','INISTP_PEAKMJD')
        else:
            pass

            
        with open('tmp.nml','w') as fout:
            print(nmltext,file=fout)
            
        os.system(f'snlc_fit.exe tmp.nml')

        return outfile

    def main(self):

        with open('opticalnir_singleband_fitresults_lowz.txt','w') as fout:
            print('# SNID z tmax tmaxerr st sterr AV AVerr Bmu Bmuerr gmu gmuerr rmu rmuerr imu imuerr Ymu Ymuerr Jmu Jmuerr Hmu Hmuerr',file=fout)

            fr = txtobj('output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',fitresheader=True)
            for j,i in enumerate(fr.CID):
                if i not in self.goodcids: continue

                outline_entries = [i,f"{fr.zHD[j]:.5f}",f"{fr.PKMJD[j]:.2f}",f"{fr.PKMJDERR[j]:.2f}",
                                   f"{fr.STRETCH[j]:.4f}",f"{fr.STRETCHERR[j]:.4f}",f"{fr.AV[j]:.4f}",f"{fr.AVERR[j]:.4f}"]
                for band in ['B','g','r','i','Yy','Jj','H']:
                    outfile = self.run_snana(
                    [i],'CSPDR3_RAISIN','kcor_CSPDR3_BD17.fits',band,
                        f"dispmodelfits/CSPDR3_RAISIN_{band}_single_{i}",clobber=False,
                        inival_av=fr.AV[j],inival_st=fr.STRETCH[j],inival_pkmjd=fr.PKMJD[j])
                    frs = txtobj(outfile+'.FITRES.TEXT',fitresheader=True)
                    if 'CID' in frs.__dict__.keys():
                        outline_entries += [f"{frs.DLMAG[0]:.5f}",f"{frs.DLMAGERR[0]:.5f}"]
                    else:
                        outline_entries += ["-99.00000","-99.00000"]

                print(" ".join(outline_entries),file=fout)

                
    def plot(self):
    
        plt.rcParams['figure.figsize'] = (7,7)
        fr = txtobj('opticalnir_singleband_fitresults_lowz.txt')
        fr.zcmb = fr.z
        fr.ebv = np.zeros(len(fr.z))
        iGood = np.array([],dtype=int)
        for i in range(len(fr.SNID)):
            fr.ebv[i] = fr.AV[i]/1.518 #lowzpars.EBVMW[lowzpars.snid == 'SN'+fr.SNID[i][2:]][0]
            if fr.z[i] > 0.035:
                iGood = np.append(iGood,i)
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]

            
        # need to shape and EBV-correct these
        for filt in filtlist:
            iGood = fr.__dict__['%smu'%filt] != -99
            zerr = 0.00083*5.0/np.log(10)*(1.0+fr.zcmb[iGood])/(fr.zcmb[iGood]*(1.0+fr.zcmb[iGood]/2.0))

            fr.__dict__['%smures'%filt] = np.array([-99.]*len(fr.__dict__['%smu'%filt]))
            fr.__dict__['%smures'%filt][iGood] = fr.__dict__['%smu'%filt][iGood] - cosmo.mu(fr.zcmb[iGood])
            fr.__dict__['%smures'%filt][iGood] -= np.average(fr.__dict__['%smures'%filt][iGood],weights=1/(fr.__dict__['%smuerr'%filt][iGood]**2.+zerr**2.))
            #for i in fr.__dict__['%smures'%filt][iGood]:
            #    if fr.__dict__['%smures'%filt] < 0:
            print(filt,'%.3f'%(weighted_avg_and_std(fr.__dict__['%smures'%filt][iGood],weights=1/(fr.__dict__['%smuerr'%filt][iGood]**2.+zerr**2.))[1]/np.std(fr.__dict__['%smures'%filt][iGood])))
            #import pdb; pdb.set_trace()
            
        covmat = np.zeros([7,7])
        for i,flt1 in enumerate(filtlist):
            for j,flt2 in enumerate(filtlist):
                imures = (fr.__dict__['%smures'%flt1] != -99) & (fr.__dict__['%smures'%flt2] != -99)
                if len(np.where(imures)[0]) > 2:
                    covmat[j,i] = np.sum(fr.__dict__['%smures'%flt1][imures]*fr.__dict__['%smures'%flt2][imures])/len(fr.__dict__['%smures'%flt1][imures])
                else:
                    covmat[j,i] = 0.0

        print(np.around(np.sqrt(np.diag(covmat)),decimals=3))
        print(np.around(cov2corr(covmat),decimals=3))
        print(nearestPD(cov2corr(covmat)))

        print(covmat)
        
        ax = plt.axes([0.1,0.07,0.8,0.7])
        cax = plt.axes([0.1,0.82,0.8,0.1])
        cb = ax.imshow(covmat)
        clb = plt.colorbar(cb,cax=cax,orientation='horizontal')
        cax.set_xlabel('$\sigma_{\mu}^2$',labelpad=-100,fontsize=20)

        ax.xaxis.set_ticks([0,1,2,3,4,5,6])
        ax.yaxis.set_ticks([0,1,2,3,4,5,6])
        ax.xaxis.set_ticklabels(['B','g','r','i','Y','J','H'])
        ax.yaxis.set_ticklabels(['B','g','r','i','Y','J','H'])
        import pdb; pdb.set_trace()

    
def covmat(samples):
    cov_shape = np.shape(samples)[1]
    chain_len = np.shape(samples)[0]
    covmat = np.zeros([cov_shape,cov_shape])
    for i in range(cov_shape):
        for j in range(cov_shape):
            covmat[j,i] = np.sum((samples[:,j]-np.mean(samples[:,j]))*(samples[:,i]-np.mean(samples[:,i])))/chain_len
    return(covmat)

def cov2corr(cov, return_std=False):
    '''convert covariance matrix to correlation matrix

    Parameters
    ----------
    cov : array_like, 2d
        covariance matrix, see Notes

    Returns
    -------
    corr : ndarray (subclass)
        correlation matrix
    return_std : bool
        If this is true then the standard deviation is also returned.
        By default only the correlation matrix is returned.

    Notes
    -----
    This function does not convert subclasses of ndarrays. This requires
    that division is defined elementwise. np.ma.array and np.matrix are allowed.

    '''
    cov = np.asanyarray(cov)
    std_ = np.sqrt(np.diag(cov))
    corr = cov / np.outer(std_, std_)
    if return_std:
        return corr, std_
    else:
        return corr

def corr2cov(corr, std):
    '''convert correlation matrix to covariance matrix given standard deviation

    Parameters
    ----------
    corr : array_like, 2d
        correlation matrix, see Notes
    std : array_like, 1d
        standard deviation

    Returns
    -------
    cov : ndarray (subclass)
        covariance matrix

    Notes
    -----
    This function does not convert subclasses of ndarrays. This requires
    that multiplication is defined elementwise. np.ma.array are allowed, but
    not matrices.

    '''
    corr = np.asanyarray(corr)
    std_ = np.asanyarray(std)
    cov = corr * np.outer(std_, std_)
    return cov
    
if __name__ == "__main__":
    #main()
    dm = dispmodel()
    dm.main()
    dm.plot()
