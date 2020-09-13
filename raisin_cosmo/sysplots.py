#!/usr/bin/env python
# D. Jones - 9/4/20

import pylab as plt
plt.ion()
import numpy as np
from txtobj import txtobj
from scipy.stats import binned_statistic

def main():
    ax1,ax2,ax3 = plt.subplot(221),plt.subplot(222),plt.subplot(223)
    ax4 = plt.subplot(224)

    lcbase = txtobj('output/cosmo_fitres/RAISIN_stat_lcparams.txt')
    lcpv = txtobj('output/cosmo_fitres/RAISIN_pecvel_lcparams.txt')
    lcpc = txtobj('output/cosmo_fitres/RAISIN_photcal_lcparams.txt')
    lcbc = txtobj('output/cosmo_fitres/RAISIN_biascor_lcparams.txt')
    lcms = txtobj('output/cosmo_fitres/RAISIN_massdivide_lcparams.txt')

    zbins = np.append(np.linspace(0.01,0.08,5),np.linspace(0.1,0.7,5))
    pvbins = binned_statistic(lcbase.zcmb,lcpv.mb-lcbase.mb,bins=zbins,statistic='median').statistic
    ax1.plot((zbins[1:]+zbins[:-1])/2.,pvbins,'o-',color='k')
    ax1.set_title('Peculiar Velocities')

    pcbins = binned_statistic(lcbase.zcmb,lcpc.mb-lcbase.mb,bins=zbins,statistic='median').statistic
    ax2.plot((zbins[1:]+zbins[:-1])/2.,pcbins,'o-',color='k')
    ax2.set_title('Phot. Cal.')

    bcbins = binned_statistic(lcbase.zcmb,lcbc.mb-lcbase.mb,bins=zbins,statistic='median').statistic
    ax3.plot((zbins[1:]+zbins[:-1])/2.,bcbins,'o-',color='k')
    ax3.set_title('Bias Corr.')

    msbins = binned_statistic(lcbase.zcmb,lcms.mb-lcbase.mb,bins=zbins,statistic='median').statistic
    ax4.plot((zbins[1:]+zbins[:-1])/2.,msbins,'o-',color='k')
    ax4.set_title('Mass Step')

    
    import pdb; pdb.set_trace()
    
if __name__ == "__main__":
    main()
