#!/usr/bin/env python

import pylab as plt
plt.ion()
import numpy as np
import cosmo
from txtobj import txtobj
import astropy.table as at

def main():

    data = at.Table.read('output/cosmo_fitres_nir/RAISIN_all_lcparams_cosmosis.txt',format='ascii',names=('name','zcmb','zhel','dz','mb','dmb','x1','dx1','color','dcolor','3rdvar','d3rdvar','cov_m_s','cov_m_c','cov_s_c','set','ra','dec','biascor'))
    
    #fr = txtobj('output/cosmo_fitres_nir/RAISIN_all_lcparams_cosmosis.txt')
    iErr = data['dmb'] < 1
    plt.errorbar(data['zcmb'][iErr],data['mb'][iErr]+19.36-cosmo.mu(data['zcmb'][iErr]),yerr=data['dmb'][iErr],fmt='o',label='RAISIN data')
    plt.axhline(0,color='k')

    z = np.linspace(0,1,100)
    plt.plot(z,cosmo.mu(z,Om=1.0,Ode=2.0)-cosmo.mu(z),label='Om=1.0,OL=2.0')
    plt.plot(z,cosmo.mu(z,Om=1.0,Ode=1.5)-cosmo.mu(z),label='Om=1.0,OL=1.5')
    plt.plot(z,cosmo.mu(z,Om=1.0,Ode=1.0)-cosmo.mu(z),label='Om=1.0,OL=1.0')
    plt.plot(z,cosmo.mu(z,Om=1.0,Ode=0.5)-cosmo.mu(z),label='Om=1.0,OL=0.5')
    
    plt.legend()
    plt.xlabel('z')
    plt.ylabel('mures')

    import pdb; pdb.set_trace()
    
if __name__ == "__main__":
    main()
