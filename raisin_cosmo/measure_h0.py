#!/usr/bin/env python
# D. Jones - 12/27/21
"""measure H0 from cepheid distances and hubble flow SNe"""

from scipy.optimize import minimize, least_squares
from txtobj import txtobj
import numpy as np

cephdict = {'2006D':(32.9196,0.122652),
            '2012fr':(31.3780,0.0558240),
            '2015F':(31.4499,0.0642833),
            '2012ht':(31.9356,0.0344569),
            '2007sr':(31.6035,0.116142),
            '2007af':(31.7717,0.0518968),
            '2009Y':(33.0939,0.205070),
            '2007A':(34.5268,0.249845)}

class H0:
    def __init__(self):
        pass

    def main(self):
        # first we need to measure H0 from the calibrator SNe

        fr = txtobj('output/RAISIN_allsys_H0.txt')
        frnames = txtobj('output/RAISIN_H0_NIR.FITRES',fitresheader=True)
        frhf = txtobj('output/cosmo_fitres_nir/RAISIN_all_lcparams_cosmosis.txt')
        
        sndist,sndisterr,cephdist,cephdisterr = np.array([]),np.array([]),np.array([]),np.array([])
        for i,c in enumerate(cephdict.keys()):
            if c not in frnames.CID:
                import pdb; pdb.set_trace()
            sndist = np.append(sndist,fr.mb[frnames.CID == c]+19.36)
            sndisterr = np.append(sndisterr,fr.dmb[frnames.CID == c])
            #sndist = np.append(sndist,frnames.DLMAG[frnames.CID == c]-19.25)
            #sndisterr = np.append(sndisterr,frnames.DLMAGERR[frnames.CID == c])
            cephdist = np.append(cephdist,cephdict[c][0])
            cephdisterr = np.append(cephdisterr,cephdict[c][1])

        def best_snmag_chi2(x):

            return -np.sum(-(sndist-cephdist-x[0])**2./(2*(sndisterr**2.+cephdisterr**2.))+\
                np.log(1/(np.sqrt(2*np.pi)*(sndisterr**2.+cephdisterr**2.))))

        def best_snmag_lsq(x):

            return (sndist-cephdist-x[0])/np.sqrt(sndisterr**2.+cephdisterr**2.)

        
        md = minimize(best_snmag_chi2,(0.0))
        md = least_squares(best_snmag_lsq,(0.0))
        J = md.jac
        cov = np.linalg.inv(J.T.dot(J))
        sn_absmag = md.x[0]
        sn_absmag_err = np.sqrt(cov[0][0])
        #import pdb; pdb.set_trace()
        zmin = 0.01
        zmax = 1.0
        sndist_full = frhf.mb[(frhf.zcmb > zmin) & (frhf.zcmb < zmax)]+19.36
        sndisterr_full = frhf.dmb[(frhf.zcmb > zmin) & (frhf.zcmb < zmax)]
        snz_full = frhf.zcmb[(frhf.zcmb > zmin) & (frhf.zcmb < zmax)]
        #sndist_full = frnames.DLMAG[frnames.zHD > zmin]-19.25
        #sndisterr_full = frnames.DLMAGERR[frnames.zHD > zmin]
        #snz_full = frnames.zHD[frnames.zHD > zmin]
        
        def best_h0_chi2(x):
            h0 = x[0]
            aB =  5*np.log10(3e5*snz_full*(1+1/2.*(1-0.55)*snz_full - 1/6.*(1-0.55-3*0.55**2.+1)*snz_full**2.))
            #- 5*np.log10(x[0])+sn_absmag
            chi2 = np.sum((sndist_full-(aB-5*np.log10(h0)+sn_absmag+25))**2./(sn_absmag_err**2.+sndisterr_full**2.))
            #import pdb; pdb.set_trace()
            return chi2

        def best_h0_lsq(x):
            logh0 = x[0]
            aB =  5*np.log10(3e5*snz_full*(1+1/2.*(1+0.55)*snz_full - 1/6.*(1+0.55-3*0.55**2.+1)*snz_full**2.))
            #- 5*np.log10(x[0])+sn_absmag
            #import pdb; pdb.set_trace()
            return (sndist_full-(aB-5*logh0+sn_absmag+25))/sndisterr_full

        def best_h0_lsq_test(x):
            logh0 = x[0]

            aB =  5*np.log10(3e5*snz_full*(1+1/2.*(1+0.55)*snz_full - 1/6.*(1-0.55-3*0.55**2.+1)*snz_full**2.))
            #- 5*np.log10(x[0])+sn_absmag
            return (sndist_full-(aB-5*logh0+sn_absmag+25))/sndisterr_full


        #aB =  np.log10(3e5*snz_full*(1+1/2.*(1+0.55)*snz_full - 1/6.*(1-0.55-3*0.55**2.+1)*snz_full**2.)) - 0.2*(sndist_full)
        
        mdh = minimize(best_h0_chi2,(70.0))
        mdh = least_squares(best_h0_lsq,(np.log10(70)))
        h0 = 10**mdh.x[0]
        J = mdh.jac
        H0magerr = np.sqrt(np.linalg.inv(J.T.dot(J))[0][0])
        h0err = np.log(10)*0.2*h0*np.sqrt(sn_absmag_err**2.+H0magerr**2.)
        print(h0,h0err)
        import pdb; pdb.set_trace()
        
        

    

if __name__ == "__main__":
    h0meas = H0()
    h0meas.main()
