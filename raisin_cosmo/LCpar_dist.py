import pandas as pd
import numpy as np
import sys,os
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.optimize import leastsq
from scipy.optimize import minimize
import argparse
import lmfit
from lmfit import report_fit
import emcee
from astropy.io import ascii
import os
import distutils.util
from lmfit import report_fit

# python raisin_cosmo/LCpar_dist.py -d output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT -s output/simdump/CSP_RAISIN_SIM_FLATDIST.DUMP -f output/fit_all/CSP_RAISIN_OPTNIR_SIM_FLATDIST/CSP_RAISIN_SIM_FLATDIST/FITOPT000.FITRES.gz --snoopy

# need combined PS1/DES sims
# snlc_sim.exe sim/flatdist/DES/sim_DES_SNOOPY.input
# snlc_sim.exe sim/flatdist/PS1/SIMGEN_PS1SPEC.INPUT
# split_and_fit.pl $RAISIN_ROOT/cosmo/fit/flatdist/DES_RAISIN_optnir.nml
# split_and_fit.pl $RAISIN_ROOT/cosmo/fit/flatdist/PS1_RAISIN_optnir.nml
# cat $SNDATA_ROOT/SIM/DES_SIM_FLATDIST_SNOOPY/DES_SIM_FLATDIST_SNOOPY.DUMP > output/simdump/HIGHZ_RAISIN_SIM_FLATDIST.DUMP
# cat $SNDATA_ROOT/SIM/PS1_SIM_FLATDIST_SNOOPY/PS1_SIM_FLATDIST_SNOOPY.DUMP >> output/simdump/HIGHZ_RAISIN_SIM_FLATDIST.DUMP
# cat output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT > output/fit_optical/HIGHZ_RAISIN_optnir.FITRES.TEXT
# cat output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT >> output/fit_optical/HIGHZ_RAISIN_optnir.FITRES.TEXT
# cat output/fit_all/DES_SIM_FLATDIST_SNOOPY/DES_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES > output/fit_all/HIGHZ_RAISIN_FLATDIST.FITRES
# cat output/fit_all/PS1_SIM_FLATDIST_SNOOPY/PS1_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES >> output/fit_all/HIGHZ_RAISIN_FLATDIST.FITRES
# go into output/simdump/HIGHZ_RAISIN_SIM_FLATDIST.DUMP,  output/fit_optical/HIGHZ_RAISIN_optnir.FITRES.TEXT, output/fit_all/HIGHZ_RAISIN_FLATDIST.FITRES, delete headers

############################
# now trying with NIR-only!!
############################
# python raisin_cosmo/LCpar_dist.py -d output/fit_nir/CSP_RAISIN_NIR_SHAPE.FITRES.TEXT -s output/simdump/CSP_RAISIN_SIM_FLATDIST.DUMP -f output/fit_all/CSP_RAISIN_NIR_SIM_FLATDIST_SHAPE/CSP_RAISIN_SIM_FLATDIST/FITOPT000.FITRES.gz -o output/fit_all/CSP_RAISIN_OPTNIR_SIM_FLATDIST/CSP_RAISIN_SIM_FLATDIST/FITOPT000.FITRES.gz --snoopy
#
# combined PS1/DES sims
# steps above and then:
# cat output/fit_nir/PS1_RAISIN_NIR_SHAPE.FITRES.TEXT > output/fit_nir/HIGHZ_RAISIN_NIR_SHAPE.FITRES.TEXT
# cat output/fit_nir/DES_RAISIN_NIR_SHAPE.FITRES.TEXT >> output/fit_nir/HIGHZ_RAISIN_NIR_SHAPE.FITRES.TEXT
# gunzip output/fit_all/PS1_SIM_NIR_FLATDIST_SNOOPY_SHAPE/PS1_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES.gz
# gunzip output/fit_all/DES_SIM_NIR_FLATDIST_SNOOPY_SHAPE/DES_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES.gz
# cat output/fit_all/PS1_SIM_NIR_FLATDIST_SNOOPY_SHAPE/PS1_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES > output/fit_nir/HIGHZ_RAISIN_FLATDIST_SIM.FITRES
# cat output/fit_all/DES_SIM_NIR_FLATDIST_SNOOPY_SHAPE/DES_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES >> output/fit_nir/HIGHZ_RAISIN_FLATDIST_SIM.FITRES
# cat output/fit_all/DES_SIM_FLATDIST_SNOOPY/DES_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES > output/fit_all/HIGHZ_RAISIN_FLATDIST.FITRES
# cat output/fit_all/PS1_SIM_FLATDIST_SNOOPY/PS1_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES >> output/fit_all/HIGHZ_RAISIN_FLATDIST.FITRES
# go into output/fit_nir/HIGHZ_RAISIN_NIR_SHAPE.FITRES.TEXT, output/fit_nir/HIGHZ_RAISIN_FLATDIST_SIM.FITRES, delete headers
# python raisin_cosmo/LCpar_dist.py -d output/fit_nir/HIGHZ_RAISIN_NIR_SHAPE.FITRES.TEXT -s output/simdump/HIGHZ_RAISIN_SIM_FLATDIST.DUMP -f output/fit_nir/HIGHZ_RAISIN_FLATDIST_SIM.FITRES -o output/fit_all/HIGHZ_RAISIN_FLATDIST.FITRES --snoopy

# this is dumb because the simulated vs. recovered populations are very different fractions
# should just take the average of both run individually
# cp $SNDATA_ROOT/SIM/DES_SIM_FLATDIST_SNOOPY/DES_SIM_FLATDIST_SNOOPY.DUMP output/simdump/DES_SIM_FLATDIST_SNOOPY.DUMP
# cp $SNDATA_ROOT/SIM/PS1_SIM_FLATDIST_SNOOPY/PS1_SIM_FLATDIST_SNOOPY.DUMP output/simdump/PS1_SIM_FLATDIST_SNOOPY.DUMP
# PS1:
# python raisin_cosmo/LCpar_dist.py -d output/fit_nir/PS1_RAISIN_NIR_SHAPE.FITRES.TEXT -s output/simdump/PS1_SIM_FLATDIST_SNOOPY.DUMP -f output/fit_all/PS1_SIM_NIR_FLATDIST_SNOOPY_SHAPE/PS1_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES -o output/fit_all/PS1_SIM_FLATDIST_SNOOPY/PS1_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES --snoopy
# DES:
# python raisin_cosmo/LCpar_dist.py -d output/fit_nir/DES_RAISIN_NIR_SHAPE.FITRES.TEXT -s output/simdump/DES_SIM_FLATDIST_SNOOPY.DUMP -f output/fit_all/DES_SIM_NIR_FLATDIST_SNOOPY_SHAPE/DES_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES -o output/fit_all/DES_SIM_FLATDIST_SNOOPY/DES_SIM_FLATDIST_SNOOPY/FITOPT000.FITRES --snoopy

#CSP
#GENPEAK_STRETCH:   0.904
#GENRANGE_STRETCH:  0.7  1.3
#GENSIGMA_STRETCH:  0.204  0.221

#STSYS
#GENPEAK_STRETCH:   0.997
#GENRANGE_STRETCH:  0.7  1.3
#GENSIGMA_STRETCH:  0.274  0.290

# PS1
#GENPEAK_STRETCH:   1.048
#GENRANGE_STRETCH:  0.7  1.3
#GENSIGMA_STRETCH:  0.197  0.212

#STSYS
#GENPEAK_STRETCH:   1.165
#GENRANGE_STRETCH:  0.7  1.3
#GENSIGMA_STRETCH:  0.268  0.278

#DES
#GENPEAK_STRETCH:   0.884
#GENRANGE_STRETCH:  0.7  1.3
#GENSIGMA_STRETCH:  0.161  0.125

#STSYS
#GENPEAK_STRETCH:   0.956
#GENRANGE_STRETCH:  0.7  1.3
#GENSIGMA_STRETCH:  0.244  0.188

# high-z avg
#NORMAL:
#GENPEAK_STRETCH:   0.966
#GENRANGE_STRETCH:  0.7  1.3
#GENSIGMA_STRETCH:  0.179  0.168
#
#STSYS
#GENPEAK_STRETCH:   1.060
#GENRANGE_STRETCH:  0.7  1.3
#GENSIGMA_STRETCH:  0.256  0.233

parser=argparse.ArgumentParser()
parser.add_argument('-d','--data_fitres', help='data survey name')
parser.add_argument('-s','--sim_dump', help='sim DUMP file w/ flag x1/c distributions')
parser.add_argument('-f','--sim_fitres', help='sim survey FITRES file')
parser.add_argument('-o','--sim_optical_fitres', help='sim survey optical FITRES file (for making cuts)')
parser.add_argument('-p','--private_data_path', default="$SNDATA_ROOT",
                    help='sim survey name w/ flag x1/c distributions')
parser.add_argument('--param', help='x1,c,STRETCH,AV')
parser.add_argument('--i_c', help='Initial guess for c', default=[0,.1,.1], nargs = '+')
parser.add_argument('--i_x1', help='Initial guess for x1', default=[0,1,1], nargs = '+')
parser.add_argument('--i_st', help='Initial guess for STRETCH', default=[1.0,.1,.1], nargs = '+')
parser.add_argument('--i_av', help='Initial guess for AV', default=[0,0,0.2], nargs = '+')
parser.add_argument('--choice', help='Which bin you wish to investigate')
parser.add_argument(
    '--snoopy', help='AV/STRETCH, not x1/c', default=False, action="store_true")
parser.add_argument(
    '--redshift', help='If you want to investigate redshift bins rather than mass', default=False)
parser.add_argument(
    '--test', help='Tests an iteration, more print statements', default=False)
parser.add_argument(
    '--cycle', help='Cycle through all the surveys', default=False)

args = parser.parse_args()
param = args.param
test = args.test
i_c = args.i_c
i_c = [float(i) for i in i_c]
i_x1 = args.i_x1
i_x1 = [float(i) for i in i_x1]
i_st = args.i_st
i_st = [float(i) for i in i_st]
i_av = args.i_av
i_av = [float(i) for i in i_av]
cycle = args.cycle

opterino = {}
opterino['maxiter'] = 1000


filename1 = os.path.expandvars(args.sim_dump)
filename2 = os.path.expandvars(args.data_fitres)
filename3 = os.path.expandvars(args.sim_fitres)
filename4 = os.path.expandvars(args.sim_optical_fitres)

def agauss(x,a,x0,sigma_l,sigma_r):
    if x < x0:
        return a*np.exp(-(x-x0)**2/(2*sigma_l**2))
    else:
        return a*np.exp(-(x-x0)**2/(2*sigma_r**2))

def tau(x,a,tau):
    return a*np.exp(-x/tau)
    
dfpre = ascii.read(filename1).to_pandas()
dfdata = ascii.read(filename2).to_pandas()
#if filename2 == 'output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT':
goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)))
#goodcids = goodcids[(goodcids != '2009al') & (goodcids != '2005el') & (goodcids != '2007as')]

debug = False
if debug:
    import glob
    import snana
    snidlist = []
    min_phase_list = []
    med_phase_list = []
    nfilt_list = []
    files = glob.glob('data/Photometry/PS1_RAISIN/*.snana.dat')
    for f in files:
        sn = snana.SuperNova(f)
        min_phase_list += [(np.min(sn.MJD[(sn.FLT == 'J') | (sn.FLT == 'H')])-sn.PEAKMJD)/(1+float(sn.REDSHIFT_HELIO.split()[0]))]
        med_phase_list += [(np.median(sn.MJD[(sn.FLT == 'J') | (sn.FLT == 'H')])-sn.PEAKMJD)/(1+float(sn.REDSHIFT_HELIO.split()[0]))]
        nfilt_list += [len(np.unique(sn.FLT[(sn.FLT == 'J') | (sn.FLT == 'H')]))]
        snidlist += [sn.SNID]
        
    files2 = glob.glob('data/Photometry/DES_RAISIN/*.snana.dat')
    for f in files2:
        sn = snana.SuperNova(f)
        min_phase_list += [(np.min(sn.MJD[(sn.FLT == 'J') | (sn.FLT == 'H')])-sn.PEAKMJD)/(1+float(sn.REDSHIFT_HELIO.split()[0]))]
        med_phase_list += [(np.median(sn.MJD[(sn.FLT == 'J') | (sn.FLT == 'H')])-sn.PEAKMJD)/(1+float(sn.REDSHIFT_HELIO.split()[0]))]
        nfilt_list += [len(np.unique(sn.FLT[(sn.FLT == 'J') | (sn.FLT == 'H')]))]
        snidlist += [sn.SNID]
        
    min_phase_list = np.array(min_phase_list)
    med_phase_list = np.array(med_phase_list)
    nfilt_list = np.array(nfilt_list)
    snidlist = np.array(snidlist)

    goodidx = np.array([],dtype=int)
    for j,i in enumerate(goodcids):
        if 'DES' not in i and 'PS' not in i: continue
        if min_phase_list[snidlist == i][0] > 10:
            goodidx = np.append(goodidx,j)
    goodcids = goodcids[goodidx]

    
iGood = np.array([],dtype=bool)
for i in range(len(dfdata)):
    if dfdata['CID'][i] in goodcids: iGood = np.append(iGood,True)
    else: iGood = np.append(iGood,False)

dfdata = dfdata[iGood]

dfpost = ascii.read(filename3).to_pandas()
dfpostopt = ascii.read(filename4).to_pandas()
idx,idxopt = np.array([],dtype=int),np.array([],dtype=int)
for j,i in enumerate(dfpostopt.CID.tolist()):
    if i in dfpost.CID.tolist():
        idx = np.append(idx,np.where((np.array(dfpost.CID) == i) & (dfpost.IDSURVEY == dfpostopt.IDSURVEY[j]))[0])
        idxopt = np.append(idxopt,j)

dfpost = dfpost.iloc[idx] #dfpost[idx]
dfpostopt = dfpostopt.iloc[idxopt] #dfpostopt[idxopt]
#import pdb; pdb.set_trace()

dfpre.CID = pd.to_numeric(dfpre.CID, errors='coerce')
dfpre = dfpre.loc[dfpre.CID == dfpre.CID]
dfpre.CID = dfpre.CID.astype(int)
#if filename2 == 'output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT':
#dfpost = dfpost[(dfpost['AV'] < 0.3*dfpost['RV']) & (dfpost['STRETCH'] > 0.75) & (dfpost['STRETCH'] < 1.185) & (dfpost['STRETCHERR'] < 1.3)]
dfpost = dfpost[(np.array(dfpostopt['AV']) < 0.3*np.array(dfpostopt['RV'])) & (np.array(dfpostopt['STRETCH']) > 0.75) & (np.array(dfpostopt['STRETCH']) < 1.185) & (np.array(dfpostopt['STRETCHERR']) < 1.3)]

# lot of nasty outliers at ~1.17-1.18 when fitting NIR-only, prevents this method from working.  Let's focus on the sample that can be accurately measured for inferring intrinsic populations
dfpost = dfpost[(dfpost['STRETCH'] < 1.165) | (dfpost['STRETCH'] > 1.185)]
#dfpost = dfpost[dfpost['FITPROB'] > 0.001]
if not args.snoopy:
    dfpre.S2c = pd.to_numeric(dfpre.S2c, errors='coerce')
    dfpre.S2x1 = pd.to_numeric(dfpre.S2x1, errors='coerce')
else:
    dfpre.STRETCH = pd.to_numeric(dfpre.STRETCH, errors='coerce')
    dfpre.AV = pd.to_numeric(dfpre.AV, errors='coerce')
#import pdb; pdb.set_trace()

def Matrix_AV_init():

    send = {}
    length = {}
    bins = np.arange(-1.0,1.0,.05)
    for i in bins:
        if i < 0.94:
            located = dfpre.loc[(dfpre.AV > i ) & (dfpre.AV <= i+0.01)]
            #print(located)
            passed = dfpost.loc[dfpost.CID.isin(np.unique(located.CID.values))]
            dfo = np.histogram(passed.AV.values, bins=bins)[0]
            send[str(np.round(i+.005,3))] = np.array(dfo)
            length[str(np.round(i+.005,3))] = len(located)
    
    for num,i in enumerate(send):
        if num == 0:
            if np.sum(length[i]) != 0:
                 combined = (send[i])/length[i]
            else: 
                 combined = (send[i])
        
        else:
            if np.sum(length[i]) != 0:
                 send[i] = (send[i])/length[i]
            combined = np.vstack((combined, send[i]))
    cI_m = np.copy(combined)
    return cI_m

def Matrix_AV(params, dfk, AV_m):

    AV_tau = params[:]
    #AV_mean, AV_l, AV_r = params
    bins = np.arange(-1.0,1.0,.05)
        
    AVdatI = np.histogram(dfk.AV.values, bins=bins)[0]

    AV = [1, AV_tau] #AV_mean, AV_l, AV_r]
    input_AV = []
    for x in ((bins[1:] + bins[:-1])/2):
        input_AV.append(tau(x, *AV))
    input_AV = np.array(input_AV)

    MP = np.matmul(input_AV.reshape([1,len(bins)-1]), AV_m)
    #import pylab as plt
    #plt.ion()
    #plt.clf()
    #plt.plot(bins[:-1],AVdatI,drawstyle='steps')
    #plt.plot(bins[:-1],input_AV*np.sum(AVdatI)/np.sum(MP),drawstyle='steps')
    #import pdb; pdb.set_trace()

    MP = MP*((np.sum(AVdatI))/np.sum(MP))
    
    Delta_AV = AVdatI - MP
    Delta_AVerr = np.copy(AVdatI)
    Delta_AVerr[Delta_AVerr == 0 ] = 1
    Delta_AVerr = np.sqrt(np.abs(Delta_AVerr))
    chi2 = []
    for i in range(len(Delta_AV)):
        temp = ((Delta_AV[i])**2)/((Delta_AVerr[i])**2)
        chi2.append(temp)
    chi2 = np.array(chi2)
    LL = -np.sum(chi2)/2.
    if LL != LL:
        LL = -np.inf
    if AV_tau > 0.6 or AV_tau < 0.03:
        LL = -np.inf
    #if (AV_l < 0.02) or (AV_r < 0.02):
    #   LL = -np.inf
    #if (AV_l > .5) or (AV_r > .5):
    #   LL = -np.inf
    #if (np.abs(AV_mean) > .3):
    #   LL = -np.inf
    return LL


def Matrix_c_init():
    #cI_mean, cI_l, cI_r = params
    send = {}
    length = {}
    bins = np.arange(-.4,.41,.02)
    for i in bins:
        if i < .4:
            located = dfpre.loc[(dfpre.S2c > i ) & (dfpre.S2c <= i+0.01)]
            #print(located)
            passed = dfpost.loc[dfpost.CID.isin(np.unique(located.CID.values))]
            dfo = np.histogram(passed.c.values, bins=bins)[0]
            send[str(np.round(i+.005,3))] = np.array(dfo)
            length[str(np.round(i+.005,3))] = len(located)
    
    for num,i in enumerate(send):
        if num == 0:
            if np.sum(length[i]) != 0:
                 combined = (send[i])/length[i]
            else: 
                 combined = (send[i])
        
        else:
            if np.sum(length[i]) != 0:
                 send[i] = (send[i])/length[i]
            combined = np.vstack((combined, send[i]))
    cI_m = np.copy(combined)
    return cI_m

def Matrix_c(params, dfk, cI_m):

    cI_mean, cI_l, cI_r = params
    bins = np.arange(-.4,.41,.02)
        
    cdatI = np.histogram(dfk.c.values, bins=bins)[0]

    cI = [1, cI_mean, cI_l, cI_r]
    input_cI = []
    for x in ((bins[1:] + bins[:-1])/2):
        input_cI.append(agauss(x, *cI))
    input_cI = np.array(input_cI)
    
    MP = np.matmul(input_cI, cI_m)
    MP = MP*((np.sum(cdatI))/np.sum(MP))
    
    Delta_c = cdatI - MP
    Delta_cerr = np.copy(cdatI)
    Delta_cerr[Delta_cerr == 0 ] = 1
    Delta_cerr = np.sqrt(np.abs(Delta_cerr))
    chi2 = []
    for i in range(len(Delta_c)):
        temp = ((Delta_c[i])**2)/((Delta_cerr[i])**2)
        chi2.append(temp)
    chi2 = np.array(chi2)
    LL = -np.sum(chi2)/2.
    if LL != LL:
        LL = -np.inf
    if (cI_l < 0) or (cI_r < 0):
        LL = -np.inf
    if (cI_l > .3) or (cI_r > .3):
        LL = -np.inf
    if (np.abs(cI_mean) > .3):
        LL = -np.inf
    return LL



def Matrix_x_init():
    send = {}
    length = {}
    bins = np.arange(-4,4.1,.2)
    for i in bins:
        if i < 4:
            located = dfpre.loc[(dfpre.S2x1 > i ) & (dfpre.S2x1 <= i+0.01)]
                            
            passed = dfpost.loc[dfpost.CID.isin(np.unique(located.CID.values))]
            dfo = np.histogram(passed.x1.values, bins=bins)[0]
            send[str(np.round(i+.05,3))] = np.array(dfo)
            length[str(np.round(i+.05,3))] = len(located)

    for num,i in enumerate(send):
        if num == 0:
            if np.sum(length[i]) != 0:
                 combined = (send[i])/length[i]
            else:
                 combined = (send[i])

        else:
            if np.sum(length[i]) != 0:
                 send[i] = (send[i])/length[i]
            combined = np.vstack((combined, send[i]))
    xI_m = np.copy(combined)
    return xI_m
    
def Matrix_x(params, dfk, xI_m):
    xI_mean, xI_l, xI_r = params
    bins = np.arange(-4,4.1,.2)


    
    xdatI = np.histogram(dfk.x1.values, bins=bins)[0]

    xI = [1, xI_mean, xI_l, xI_r]
    input_xI = []
    for x in ((bins[1:] + bins[:-1])/2):
        input_xI.append(agauss(x, *xI))
    input_xI = np.array(input_xI)
    
    MP = np.matmul(input_xI, xI_m)
    MP = MP*((np.sum(xdatI))/np.sum(MP))
    
    Delta_x = xdatI - MP
    Delta_xerr = np.copy(xdatI)
    Delta_xerr[Delta_xerr == 0 ] = 1
    Delta_xerr = np.sqrt(np.abs(Delta_xerr))
    chi2 = []
    for i in range(len(Delta_x)):
        temp = ((Delta_x[i])**2)/((Delta_xerr[i])**2)
        chi2.append(temp)
    chi2 = np.array(chi2)
    LL = -np.sum(chi2)/2.
    if LL != LL:
        LL = -np.inf
    if (xI_l < 0) or (xI_r < 0):
        LL = -np.inf
    if (xI_l > 3) or (xI_r > 3):
        LL = -np.inf
    if (np.abs(xI_mean) > 3):
        LL = -np.inf
    return LL

def Matrix_stretch_init():
    send = {}
    length = {}
    bins = np.arange(0.6,1.4,.05)
    for i in bins:
        if i < np.max(bins):
            located = dfpre.loc[(dfpre.STRETCH > i ) & (dfpre.STRETCH <= i+0.01)]
                            
            passed = dfpost.loc[dfpost.CID.isin(np.unique(located.CID.values))]
            dfo = np.histogram(passed.STRETCH.values, bins=bins)[0]
            send[str(np.round(i+.05,3))] = np.array(dfo)
            length[str(np.round(i+.05,3))] = len(located)
            #import pdb; pdb.set_trace()
    for num,i in enumerate(send):
        if num == 0:
            if np.sum(length[i]) != 0:
                 combined = (send[i])/length[i]
            else:
                 combined = (send[i])

        else:
            if np.sum(length[i]) != 0:
                 send[i] = (send[i])/length[i]
            combined = np.vstack((combined, send[i]))
    xI_m = np.copy(combined)

    return xI_m
    
def Matrix_stretch(params, dfk, xI_m):
    xI_mean, xI_l, xI_r = params
    bins = np.arange(0.6,1.4,.05)


    
    xdatI = np.histogram(dfk.STRETCH.values, bins=bins)[0]

    xI = [1, xI_mean, xI_l, xI_r]
    input_xI = []
    for x in ((bins[1:] + bins[:-1])/2):
        input_xI.append(agauss(x, *xI))
    input_xI = np.array(input_xI)
    
    MP = np.matmul(input_xI, xI_m)
    #import pylab as plt
    #plt.ion()
    #plt.clf()
    #plt.plot(bins[:-1],xdatI,drawstyle='steps')
    #plt.plot(bins[:-1],input_xI*np.sum(xdatI)/np.sum(MP),drawstyle='steps')
    #import pdb; pdb.set_trace()
    MP = MP*((np.sum(xdatI))/np.sum(MP))
    
    Delta_x = xdatI - MP
    Delta_xerr = np.copy(xdatI)
    Delta_xerr[Delta_xerr == 0 ] = 1
    Delta_xerr = np.sqrt(np.abs(Delta_xerr))
    chi2 = []
    for i in range(len(Delta_x)):
        temp = ((Delta_x[i])**2)/((Delta_xerr[i])**2)
        chi2.append(temp)
    chi2 = np.array(chi2)
    LL = -np.sum(chi2)/2.
    #print(LL)
    #if LL != LL:
    #import pdb; pdb.set_trace()
    if LL != LL:
        LL = -np.inf
    if (xI_l < 0.02) or (xI_r < 0.02):
        LL = -np.inf
    if (xI_l > 0.3) or (xI_r > 0.3):
        LL = -np.inf
    if (np.abs(xI_mean) > 1.185):
        LL = -np.inf
    if (np.abs(xI_mean) < 0.75):
        LL = -np.inf

    return LL


def Optimiser():
    cI_m = Matrix_c_init()
    nwalkers = 7
    ndim = 3
    p0 = np.random.rand(nwalkers, ndim)
    p0 = p0/100
    p0 = np.abs(p0)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix_c, args=[dfk, cI_m])
    state = sampler.run_mcmc(p0, 1000, progress=True)
    sampler.reset()
    sampler.run_mcmc(state, 10000, progress=True)
    print(np.mean(sampler.acceptance_fraction))
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    labels = ["mean", "stdl", "stdr"]
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number");
    plt.savefig('/home/bap37/Documents/Cosmology/MASS/TEST.pdf', format='pdf')
    
#Now we need to cycle through masses
class Optimiser_Mass:

    def main(self):
        stdls = []
        track = []
        stdrs = []
        means = []
        mass = []
        c = 10
        massstepper = np.arange(6,14,.2)
        m=10
        for param in ['AV','STRETCH']:
            fbf=False
            dfk = dfdata
            print(len(dfk))
            min_count = 15
            if (len(dfk) > min_count) & fbf==False:
                track.append(m)
                fbf = True
            if len(dfk) < min_count:
                pass
            else:
                if param == 'STRETCH':
                    xI_m = Matrix_stretch_init()
                    nwalkers = 7
                    ndim = 3
                    p0 = np.zeros([nwalkers,ndim])
                    p0[:,0] = np.random.normal(size=nwalkers,loc=1.0,scale=0.1)
                    p0[:,1] = np.random.normal(size=nwalkers,loc=0.1,scale=0.1)
                    p0[:,2] = np.random.normal(size=nwalkers,loc=0.1,scale=0.1)
                    #p0 = np.random.rand(nwalkers, ndim)+0.5
                    #p0 = p0
                    p0 = np.abs(p0)
                    #from scipy.optimize import minimize
                    #nll = lambda *args: -Matrix_stretch(*args)
                    #md = minimize(Matrix_stretch,(1,0.1,0.1),args=(dfk,xI_m),bounds=((0.8,1.2),(0.05,0.3),(0.05,0.3)))
                    sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix_stretch, args=[dfk, xI_m])
                    state = sampler.run_mcmc(p0, 2000, progress=True)
                    samples = sampler.get_chain()
                    sampler.reset()
                    sampler.run_mcmc(state, 10000, progress=True)            
                    samples = sampler.get_chain()
                    mass.append(m)
                    print(f"stretch: {np.mean(samples[:,:,0])} +/- {np.std(samples[:,:,0])}")
                    print(f"STD left: {np.mean(samples[:,:,1])} +/- {np.std(samples[:,:,1])}")
                    print(f"STD right: {np.mean(samples[:,:,2])} +/- {np.std(samples[:,:,2])}")
                    means.append([np.mean(samples[:,:,0]), np.std(samples[:,:,0])])
                    stdls.append([np.mean(samples[:,:,1]), np.std(samples[:,:,1])])
                    stdrs.append([np.mean(samples[:,:,2]), np.std(samples[:,:,2])])
                    
                elif param == 'AV':
                    cI_m = Matrix_AV_init()
                    nwalkers = 7
                    ndim = 1
                    p0 = np.zeros([nwalkers,ndim])
                    p0[:,0] = np.random.normal(size=nwalkers,loc=0.1,scale=0.1)
                    #p0 = np.random.rand(nwalkers, ndim)
                    #p0 = p0/1000
                    p0 = np.abs(p0)
                    #import pdb; pdb.set_trace()
                    sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix_AV, args=[dfk, cI_m])
                    try:
                        state = sampler.run_mcmc(p0, 2000, progress=True)
                    except ValueError:
                        state = sampler.run_mcmc(p0, 2000, progress=True)
                    print(np.mean(sampler.acceptance_fraction))
                    sampler.reset()

                    try:
                        sampler.run_mcmc(state, 10000, progress=True)
                    except ValueError:
                        sampler.run_mcmc(state, 10000, progress=True)
                    samples = sampler.get_chain()
                    mass.append(m)
                    print(f"AV tau: {np.mean(samples[:,:,0])} +/- {np.std(samples[:,:,0])}")
                    means.append([np.mean(samples[:,:,0]), np.std(samples[:,:,0])])
                    #stdls.append([np.mean(samples[:,:,1]), np.std(samples[:,:,1])])
                    #stdrs.append([np.mean(samples[:,:,2]), np.std(samples[:,:,2])])
                    
        print(means)
        print(stdls)
        print(stdrs)


        print(f"""NORMAL:
GENPEAK_STRETCH:   {means[1][0]:.3f}
GENRANGE_STRETCH:  0.7  1.3
GENSIGMA_STRETCH:  {stdls[0][0]:.3f}  {stdrs[0][0]:.3f}

GENPEAK_RV:         1.518             # mean RV to generate
GENSIGMA_RV:        0.001  0.001      # lower,upper guassian-sigmas
GENRANGE_AV:        0.0   2.00        # CCM89 extinc param range
GENTAU_AV:          {means[0][0]:.3f}            # dN/dAV = exp(-AV/tau)

STSYS
GENPEAK_STRETCH:   {means[1][0]+means[1][1]:.3f}
GENRANGE_STRETCH:  0.7  1.3
GENSIGMA_STRETCH:  {stdls[0][0]+stdls[0][1]:.3f}  {stdrs[0][0]+stdrs[0][1]:.3f}

GENPEAK_RV:         1.518             # mean RV to generate
GENSIGMA_RV:        0.001  0.001      # lower,upper guassian-sigmas
GENRANGE_AV:        0.0   2.00        # CCM89 extinc param range
GENTAU_AV:          {means[0][0]:.3f}            # dN/dAV = exp(-AV/tau)

AVSYS:
GENPEAK_STRETCH:   {means[1][0]:.3f}
GENRANGE_STRETCH:  0.7  1.3
GENSIGMA_STRETCH:  {stdls[0][0]:.3f}  {stdrs[0][0]:.3f}

GENPEAK_RV:         1.518             # mean RV to generate
GENSIGMA_RV:        0.001  0.001      # lower,upper guassian-sigmas
GENRANGE_AV:        0.0   2.00        # CCM89 extinc param range
GENTAU_AV:          {means[0][0]+means[0][1]:.3f}            # dN/dAV = exp(-AV/tau)


""")
        
        sys.exit()
        data_survey = self.options.data_survey
        with open(f'{self.options.data_survey}_{param}_masses.txt',"w+") as fout:
            for m in mass:
                fout.write(str(m)+'\n')
        with open(f'{self.options.data_survey}_{param}_means.txt',"w+") as fout:
            for f in means:
                fout.write(str(f)+'\n')
        with open(f'{self.options.data_survey}_{param}_stdls.txt',"w+") as fout:
            for f in stdls:
                fout.write(str(f)+'\n')
        with open(f'{self.options.data_survey}_{param}_stdrs.txt',"w+") as fout:
            for f in stdls:
                fout.write(str(f)+'\n')

        return

def average_results(stmean=[1.046,0.880],ststd1=[0.194,0.158],ststd2=[0.211,0.129],stsysmean=[1.162,0.953],stsysstd1=[0.266,0.239],stsysstd2=[0.277,0.193]):
    # just need a flat average for PS1/DES

    stmean = np.average(stmean)
    ststd1 = np.average(ststd1)
    ststd2 = np.average(ststd2)
    stsysmean = np.average(stsysmean)
    stsysstd1 = np.average(stsysstd1)
    stsysstd2 = np.average(stsysstd2)
    
    print(f"""NORMAL:
GENPEAK_STRETCH:   {stmean:.3f}
GENRANGE_STRETCH:  0.7  1.3
GENSIGMA_STRETCH:  {ststd1:.3f}  {ststd2:.3f}

STSYS
GENPEAK_STRETCH:   {stsysmean:.3f}
GENRANGE_STRETCH:  0.7  1.3
GENSIGMA_STRETCH:  {stsysstd1:.3f}  {stsysstd2:.3f}

""")

    
if __name__ == "__main__": 
    import distutils.util

    #om = Optimiser_Mass()
    #om.options = args
    #om.main()
    
    average_results()

## NIR_shape failed
# DES opt failed
