#!/usr/bin/env python

from txtobj import txtobj
import numpy as np
import cosmo
import matplotlib.pyplot as plt
import numpy.linalg as la
plt.ion()
filtlist = 'BgriYJH'
lowzpars = txtobj('../snoopy_nir_pipeline/GPfit/lcparams_avelino19.txt')

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
    fr = txtobj('../snoopy_nir_pipeline/dispmodel/opticalnir_maxmodel_fitresults_lowz.txt')
    #fr = txtobj('opticalnir_maxmodel_fitresults_lowz.txt')
    fr.zcmb = np.zeros(len(fr.z))
    fr.ebv = np.zeros(len(fr.z))
    for i in range(len(fr.SNID)):
        fr.zcmb[i] = lowzpars.z_cmb[lowzpars.snid == 'SN'+fr.SNID[i][2:]][0]
        fr.ebv[i] = lowzpars.EBVMW[lowzpars.snid == 'SN'+fr.SNID[i][2:]][0]
    
    for filt in filtlist:
        iGood = fr.__dict__['%smu'%filt] != -99
        fr.__dict__['%smures'%filt] = np.array([-99.]*len(fr.__dict__['%smu'%filt]))
        fr.__dict__['%smures'%filt][iGood] = fr.__dict__['%smu'%filt][iGood] - cosmo.mu(fr.zcmb[iGood])
        fr.__dict__['%smures'%filt][iGood] -= np.average(fr.__dict__['%smures'%filt][iGood],weights=1/fr.__dict__['%smuerr'%filt][iGood]**2.)
            
    covmat = np.zeros([7,7])
    for i,flt1 in enumerate(filtlist):
        for j,flt2 in enumerate(filtlist):
            imures = (fr.__dict__['%smures'%flt1] != -99) & (fr.__dict__['%smures'%flt2] != -99)
            if len(np.where(imures)[0]) > 2:
                covmat[j,i] = np.sum(fr.__dict__['%smures'%flt1][imures]*fr.__dict__['%smures'%flt2][imures])/len(fr.__dict__['%smures'%flt1][imures])
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
    main()
