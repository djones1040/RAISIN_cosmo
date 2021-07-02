#!/usr/bin/env python
# D. Jones - 9/16/20
import numpy as np
import pylab as plt
import os
from txtobj import txtobj
import glob
plt.rcParams['figure.figsize'] =(12,12)
plt.ion()
from raisin_cosmo import ovdatamc
import palettable
#from palettable.wesanderson import Aquatic1_5 as palettable_color
from palettable.colorbrewer.qualitative import Dark2_8 as palettable_color
from scipy.stats import binned_statistic
import getmu
import copy

_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)))

_outdirs = ['$RAISIN_ROOT/cosmo/output/fit_nir_sys/CSP_RAISIN',
            '$RAISIN_ROOT/cosmo/output/fit_nir_sys/PS1_RAISIN',
            '$RAISIN_ROOT/cosmo/output/fit_nir_sys/DES_RAISIN']

_nirsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN_NIR_SIM',
                     '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN_NIR_SIM',
                     '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN_NIR_SIM']
_nirsimfitreslist = [os.path.expandvars(filepath) for filepath in _nirsimfitreslist]

_opticalnirsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_all/CSP_RAISIN_OPTNIR_SIM',
                            '$RAISIN_ROOT/cosmo/output/fit_all/PS1_RAISIN_OPTNIR_SIM',
                            '$RAISIN_ROOT/cosmo/output/fit_all/DES_RAISIN_OPTNIR_SIM']
_opticalnirsimfitreslist = [os.path.expandvars(filepath) for filepath in _opticalnirsimfitreslist]


_opticalsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_all/CSP_RAISIN_OPT_SIM',
                         '$RAISIN_ROOT/cosmo/output/fit_all/PS1_RAISIN_OPT_SIM',
                         '$RAISIN_ROOT/cosmo/output/fit_all/DES_RAISIN_OPT_SIM']
_opticalsimfitreslist = [os.path.expandvars(filepath) for filepath in _opticalsimfitreslist]


_nirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN.FITRES.TEXT',
                      '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN.FITRES.TEXT',
                      '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN.FITRES.TEXT']
_nirdatafitreslist = [os.path.expandvars(filepath) for filepath in _nirdatafitreslist]

_opticaldatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optical.FITRES.TEXT',
                          '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optical.FITRES.TEXT',
                          '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optical.FITRES.TEXT']
_opticaldatafitreslist = [os.path.expandvars(filepath) for filepath in _opticaldatafitreslist]


_opticalnirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',
                             '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT',
                             '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT']
_opticalnirdatafitreslist = [os.path.expandvars(filepath) for filepath in _opticalnirdatafitreslist]

#_g10fitreslist = ['output/fit_all/CSP_RAISIN_OPTNIR_SIM_G10/CSP_RAISIN_SIM_G10/FITOPT000.FITRES',
#                  'output/fit_all/PS1_RAISIN_OPTNIR_SIM_G10/PS1_RAISIN_SIM_G10/FITOPT000.FITRES',
#                  'output/fit_all/DES_RAISIN_OPTNIR_SIM_G10/DES_RAISIN_SIM_G10/FITOPT000.FITRES']
_g10fitreslist = ['output/fit_all/CSP_RAISIN_OPTNIR_SIM_G10/CSP_RAISIN_SIM_G10/FITOPT000.FITRES',
                  'output/fit_all/PS1_RAISIN_OPTNIR_SIM_G10/PS1_RAISIN_SIM_G10/FITOPT000.FITRES',
                  'output/fit_all/DES_RAISIN_OPTNIR_SIM_G10/DES_RAISIN_SIM_G10/FITOPT000.FITRES']

def apply_all_cuts(fr,fropt,restrict_to_good_list=False):

    # AV
    iGoodAV = np.zeros(len(fr.CID),dtype=bool)
    for j,i in enumerate(fr.CID):
        if i in fropt.CID and fropt.AV[fropt.CID == i] < 0.3*fropt.RV[fropt.CID == i]:
            iGoodAV[j] = True

    # reasonable stretch
    iGoodSt = np.zeros(len(fr.CID),dtype=bool)
    for j,i in enumerate(fr.CID):
        if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > 0.8 and fropt.STRETCH[fropt.CID == i] < 1.3: #175:
            iGoodSt[j] = True

    #iGoodSt = (fropt.STRETCH > 0.8) & (fropt.STRETCH < 1.3)

    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iGoodAV & iGoodSt]

    #restrict_to_good_list = False
    if restrict_to_good_list:
        # then case-by-case removal of bad stuff
        # because some RAISIN things have bad photometry
        iGood = np.array([],dtype=int)
        for j,i in enumerate(fr.CID):
            if i in _goodcids: iGood = np.append(iGood,j)
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]

    return fr


def main():

    plt.subplots_adjust(wspace=0.0,hspace=0.0)
    
    for nirdatadir,opticalnirdatafitres,opticaldatafitres,nirsimfitres,opticalnirsimfitres,opticalsimfitres,name,idx in zip(
            _outdirs,_opticalnirdatafitreslist,_opticaldatafitreslist,
            _nirsimfitreslist,_opticalnirsimfitreslist,_opticalsimfitreslist,['CSP','PS1','DES'],range(3)):

        if name == 'CSP':
            simname = 'CSP_RAISIN_SIM'
            goodcids = np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)
        elif name == 'PS1':
            simname = 'PS1_RAISIN_SIM'
            goodcids = np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)
        elif name == 'DES':
            simname = 'DES_RAISIN_SIM'
            goodcids = np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)
        
        ax1,ax2,ax3,ax4 = plt.subplot(3,4,4*idx+1),plt.subplot(3,4,4*idx+2),plt.subplot(3,4,4*idx+3),plt.subplot(3,4,4*idx+4)
        if name == 'CSP': plt.gcf().text(0.03,0.77,name,ha='center',va='center',rotation=90,fontsize=20)
        elif name == 'PS1': plt.gcf().text(0.03,0.53,name,ha='center',va='center',rotation=90,fontsize=20)
        elif name == 'DES': plt.gcf().text(0.03,0.28,name,ha='center',va='center',rotation=90,fontsize=20)
        
        fropt = txtobj(opticalnirdatafitres,fitresheader=True)

        print(glob.glob(os.path.expandvars(f'{opticalnirsimfitres}/{simname}/FITOPT000.FITRES'))[0])
        froptsim = txtobj(glob.glob(os.path.expandvars(f'{opticalnirsimfitres}/{simname}/FITOPT000.FITRES'))[0],
                          fitresheader=True)

        froptsimav = txtobj(glob.glob(os.path.expandvars(f'{opticalnirsimfitres}/{simname}_AVSYS/FITOPT000.FITRES'))[0],
                            fitresheader=True)
        froptsimst = txtobj(glob.glob(os.path.expandvars(f'{opticalnirsimfitres}/{simname}_STSYS/FITOPT000.FITRES'))[0],
                            fitresheader=True)

        iGood = np.array([],dtype=bool)
        for j,i in enumerate(fropt.CID):
            if i in goodcids: iGood = np.append(iGood,True)
            else: iGood = np.append(iGood,False)
        for k in fropt.__dict__.keys():
            fropt.__dict__[k] = fropt.__dict__[k][iGood]
            
        froptsim = apply_all_cuts(froptsim,froptsim,restrict_to_good_list=False)
        froptsimav = apply_all_cuts(froptsimav,froptsimav,restrict_to_good_list=False)
        froptsimst = apply_all_cuts(froptsimst,froptsimst,restrict_to_good_list=False)

        # OPT+NIR
        hist = ovdatamc.ovhist()
        parser = hist.add_options(usage='')
        options,  args = parser.parse_args()
        hist.options = options
        hist.options.journal = True
        if name != 'CSP': hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.8,1.3),('AV',-0.5,0.93),('SNRMAX1',0,60)]
        else: hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.8,1.3),('AV',-0.5,0.93),('SNRMAX1',0,150)]
        hist.options.nbins = 10
        hist.options.clobber = True
        hist.options.outfile = 'figs/sim_%s.png'%name
        hist.options.histvar = ['MURES','AV','STRETCH']#,'SNRMAX1']
        hist.main(data=fropt,sim=froptsim,axes=[ax1,ax2,ax3],letters=False,ncol=4)
        hist.main(data=fropt,sim=froptsimav,axes=[ax1,ax2,ax3],letters=False,simcolor='C0',
                  simlabel='$A_V$ 1$\sigma$',datalabel=None,simls='--',showdata=False,showmeanstd=False,ncol=4)
        hist.main(data=fropt,sim=froptsimst,axes=[ax1,ax2,ax3],letters=False,simcolor='C1',
                  simlabel='$s_{BV}$ 1$\sigma$',datalabel=None,simls='-.',showdata=False,showmeanstd=False,ncol=4)

        # OPT only
        fropt = txtobj(opticaldatafitres,fitresheader=True)

        print(glob.glob(os.path.expandvars(f'{opticalsimfitres}/{simname}/FITOPT000.FITRES'))[0])
        froptsim = txtobj(glob.glob(os.path.expandvars(f'{opticalsimfitres}/{simname}/FITOPT000.FITRES'))[0],
                          fitresheader=True)

        froptsimav = txtobj(glob.glob(os.path.expandvars(f'{opticalsimfitres}/{simname}_AVSYS/FITOPT000.FITRES'))[0],
                            fitresheader=True)
        froptsimst = txtobj(glob.glob(os.path.expandvars(f'{opticalsimfitres}/{simname}_STSYS/FITOPT000.FITRES'))[0],
                            fitresheader=True)

        iGood = np.array([],dtype=bool)
        for j,i in enumerate(fropt.CID):
            if i in goodcids: iGood = np.append(iGood,True)
            else: iGood = np.append(iGood,False)
        for k in fropt.__dict__.keys():
            fropt.__dict__[k] = fropt.__dict__[k][iGood]
            
        froptsim = apply_all_cuts(froptsim,froptsim,restrict_to_good_list=False)
        froptsimav = apply_all_cuts(froptsimav,froptsimav,restrict_to_good_list=False)
        froptsimst = apply_all_cuts(froptsimst,froptsimst,restrict_to_good_list=False)


        hist = ovdatamc.ovhist()
        parser = hist.add_options(usage='')
        options,  args = parser.parse_args()
        hist.options = options
        hist.options.journal = True
        if name != 'CSP': hist.options.cutwin = [('SNRMAX1',0,60)]
        else: hist.options.cutwin = [('SNRMAX1',0,150)]
        hist.options.nbins = 10
        hist.options.clobber = True
        hist.options.outfile = 'figs/sim_%s.png'%name
        hist.options.histvar = ['SNRMAX1']
        hist.main(data=fropt,sim=froptsim,axes=[ax4],letters=False,showlegend=False,ylabel=False)
        hist.main(data=fropt,sim=froptsimav,axes=[ax4],letters=False,simcolor='C0',
                  simlabel='$A_V$ 1$\sigma$',datalabel=None,simls='--',showdata=False,showmeanstd=False,showlegend=False,ylabel=False)
        hist.main(data=fropt,sim=froptsimst,axes=[ax4],letters=False,simcolor='C1',
                  simlabel='$s_{BV}$ 1$\sigma$',datalabel=None,simls='-.',showdata=False,showmeanstd=False,showlegend=False,ylabel=False)

        for ax in [ax1,ax2,ax3,ax4]:
            if name == 'CSP': ax.set_ylim([0,55])
            else: ax.set_ylim([0,12])
            ax.tick_params(top="off",bottom="on",direction="inout",length=8, width=2)
        ax1.set_xlim([-0.7,0.6])
        ax2.set_xlim([-0.25,0.9])
        ax3.set_xlim([0.75,1.3])
        ax4.set_xlim([0,100])
        for ax in [ax2,ax3,ax4]:
            ax.yaxis.set_ticklabels([])
        if name != 'DES':
            for ax in [ax2,ax3,ax4]:
                ax.xaxis.set_ticklabels([])
        if name != 'CSP':
            ax1.get_legend().remove()
        else:
            ax1.get_legend()._loc_real = 3
            ax1.get_legend().set_bbox_to_anchor([0.949,1.0], transform = ax1.transAxes)
            #ax1.get_legend()._ncol = 4
            #import pdb; pdb.set_trace()
        #if name == 'CSP': import pdb; pdb.set_trace()
    
    import pdb; pdb.set_trace()

def biascor():
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(223)
    ax3 = plt.subplot(224)
    for ax in [ax1,ax2,ax3]:
        ax.set_prop_cycle('color', palettable_color.mpl_colors)

    frpanlowz = txtobj('data/PanSims/Pan_LOWZ_G10.FITRES',fitresheader=True)
    frpanps1 = txtobj('data/PanSims/Pan_PS1_G10.FITRES',fitresheader=True)
    
    for frfile,froptfile,g10file,color,name in zip(
            _nirsimfitreslist,_opticalnirsimfitreslist,_g10fitreslist,palettable_color.mpl_colors,['CSP','PS1','DES']):
        if name == 'CSP': simname = 'CSP_RAISIN_SIM'
        elif name == 'PS1': simname = 'PS1_RAISIN_SIM'
        elif name == 'DES': simname = 'DES_RAISIN_SIM'

        frsim = txtobj(f'{frfile}/{simname}/FITOPT000.FITRES',fitresheader=True)
        froptsim = txtobj(f'{froptfile}/{simname}/FITOPT000.FITRES',fitresheader=True)
        frsim = apply_all_cuts(frsim,froptsim,restrict_to_good_list=False)
        froptsim = apply_all_cuts(froptsim,froptsim,restrict_to_good_list=False)

        frg = txtobj(g10file,fitresheader=True)
        frg = apply_all_cuts(frg,frg,restrict_to_good_list=False) #getmu.mkcuts(frg)
        frg.SIM_DLMAG = frg.SIM_mB + frg.SIM_alpha*frg.SIM_x1 - frg.SIM_beta*frg.SIM_c + 19.36
        iFP = frg.FITPROB > 1e-3
        for k in frg.__dict__.keys():
            frg.__dict__[k] = frg.__dict__[k][iFP]

        frpanlowz = getmu.mkcuts(frpanlowz)
        frpanps1 = getmu.mkcuts(frpanps1)
        frpanlowz = getmu.getmu(frpanlowz)
        frpanps1 = getmu.getmu(frpanps1)

        frpanlowz.DLMAG = frpanlowz.mB + frpanlowz.SIM_alpha*frpanlowz.x1 - frpanlowz.SIM_beta*frpanlowz.c + 19.36
        frpanlowz.DLMAGERR = frpanlowz.muerr
        frpanlowz.SIM_DLMAG = frpanlowz.SIM_mB + frpanlowz.SIM_alpha*frpanlowz.SIM_x1 - frpanlowz.SIM_beta*frpanlowz.SIM_c + 19.36
        frpanps1.DLMAG = frpanps1.mB + frpanps1.SIM_alpha*frpanps1.x1 - frpanps1.SIM_beta*frpanps1.c + 19.36
        frpanps1.SIM_DLMAG = frpanps1.SIM_mB + frpanps1.SIM_alpha*frpanps1.SIM_x1 - frpanps1.SIM_beta*frpanps1.SIM_c + 19.36
        frpanps1.DLMAGERR = frpanps1.muerr
        
        if name == 'CSP':
            frgtot = copy.deepcopy(frg)
            frsimtot = copy.deepcopy(frsim)
            froptsimtot = copy.deepcopy(froptsim)
        else:
            if name == 'PS1':
                froptsimdes = txtobj(f'{_nirsimfitreslist[-1]}/DES_RAISIN_SIM/FITOPT000.FITRES',fitresheader=True)
                froptsimdes = apply_all_cuts(froptsimdes,froptsim,restrict_to_good_list=False)
                frsimdes = txtobj(f'{_nirsimfitreslist[-1]}/DES_RAISIN_SIM/FITOPT000.FITRES',fitresheader=True)
                frsimdes = apply_all_cuts(frsimdes,froptsimdes,restrict_to_good_list=False)
                frgdes = txtobj(f'{_g10fitreslist[-1]}',fitresheader=True)
                frgdes = apply_all_cuts(frgdes,frgdes,restrict_to_good_list=False)
                iFP = frgdes.FITPROB > 1e-3
                for k in frgdes.__dict__.keys():
                    frgdes.__dict__[k] = frgdes.__dict__[k][iFP]
                
                from random import sample
                iRand = sample(range(len(frsim.CID)),len(frsimdes.CID))
                iRandOpt = sample(range(len(froptsim.CID)),len(froptsimdes.CID))
                iRandg = sample(range(len(frg.CID)),len(frgdes.CID))
                
                for k in frgtot.__dict__.keys():
                    frgtot.__dict__[k] = np.append(frgtot.__dict__[k],frg.__dict__[k][iRandg])
                for k in frsimtot.__dict__.keys():
                    frsimtot.__dict__[k] = np.append(frsimtot.__dict__[k],frsim.__dict__[k][iRand])
                for k in froptsimtot.__dict__.keys():
                    froptsimtot.__dict__[k] = np.append(froptsimtot.__dict__[k],froptsim.__dict__[k][iRandOpt])
                
            else:

                for k in frgtot.__dict__.keys():
                    frgtot.__dict__[k] = np.append(frgtot.__dict__[k],frg.__dict__[k])
                for k in frsimtot.__dict__.keys():
                    frsimtot.__dict__[k] = np.append(frsimtot.__dict__[k],frsim.__dict__[k])
                for k in froptsimtot.__dict__.keys():
                    froptsimtot.__dict__[k] = np.append(froptsimtot.__dict__[k],froptsim.__dict__[k])
        #import pdb; pdb.set_trace()
    #frg = getmu.getmu(frg)
    #frg.mubias = (frg.mB-frg.SIM_mB) + 0.14*(frg.x1-frg.SIM_x1) - 3.1*(frg.c - frg.SIM_c)

    def weighted_avg(x,fr=None,var='DLMAG'):
        average = np.average(fr.__dict__[var][x]-fr.__dict__['SIM_'+var][x],
                             weights=1/fr.__dict__[var+'ERR'][x]**2.)
        variance = np.average((fr.__dict__[var][x]-fr.__dict__['SIM_'+var][x]-average)**2,
                              weights=1/fr.__dict__[var+'ERR'][x]**2.)
        return average

    def weighted_std(x,fr=None,var='DLMAG'):
        average = np.average(fr.__dict__[var][x]-fr.__dict__['SIM_'+var][x],
                             weights=1/fr.__dict__[var+'ERR'][x]**2.)
        variance = np.average((fr.__dict__[var][x]-fr.__dict__['SIM_'+var][x]-average)**2,
                              weights=1/fr.__dict__[var+'ERR'][x]**2.)
        return np.sqrt(variance/float(len(x)))

    #def weighted_avg_g10(x,fr=None):
    #    average = np.average(fr.mubias[x], weights=1/fr.muerr[x]**2.)
    #    variance = np.average((fr.mubias[x]-average)**2, weights=1/fr.muerr[x]**2.)
    #    return average

    #def weighted_std_g10(x,fr=None):
    #    average = np.average(fr.mubias[x], weights=1/fr.muerr[x]**2.)
    #    variance = np.average((fr.mubias[x]-average)**2, weights=1/fr.muerr[x]**2.)
    #    return np.sqrt(variance/float(len(x)))


    zbins = np.linspace(0.01,0.8,30)
    frsimtot.AVERR[:] = 0.01
    frsimtot.STRETCHERR[:] = 0.01
    for var,ax in zip(['DLMAG','AV','STRETCH'],[ax1,ax2,ax3]):
        delmusim = binned_statistic(frsimtot.zCMB,range(len(frsimtot.zCMB)),bins=zbins,
                                    statistic=lambda values: weighted_avg(values,frsimtot,var)).statistic
        delmusimerr = binned_statistic(frsimtot.zCMB,range(len(frsimtot.zCMB)),bins=zbins,
                                       statistic=lambda values: weighted_std(values,frsimtot,var)).statistic

        delmuoptsim = binned_statistic(froptsimtot.zCMB,range(len(froptsimtot.DLMAG)),bins=zbins,
                                       statistic=lambda values: weighted_avg(values,froptsimtot,var)).statistic
        delmuoptsimerr = binned_statistic(froptsimtot.zCMB,range(len(froptsimtot.DLMAG)),bins=zbins,
                                          statistic=lambda values: weighted_std(values,froptsimtot,var)).statistic

        if var == 'DLMAG':
            delmulowz = binned_statistic(frpanlowz.zCMB,range(len(frpanlowz.zCMB)),bins=zbins,
                                         statistic=lambda values: weighted_avg(values,frpanlowz,var)).statistic
            delmulowzerr = binned_statistic(frpanlowz.zCMB,range(len(frpanlowz.zCMB)),bins=zbins,
                                            statistic=lambda values: weighted_std(values,frpanlowz,var)).statistic
            delmups1 = binned_statistic(frpanps1.zCMB,range(len(frpanps1.zCMB)),bins=zbins,
                                        statistic=lambda values: weighted_avg(values,frpanps1,var)).statistic
            delmups1err = binned_statistic(frpanps1.zCMB,range(len(frpanps1.zCMB)),bins=zbins,
                                           statistic=lambda values: weighted_std(values,frpanps1,var)).statistic

        #    delmug10sim = binned_statistic(frgtot.zCMB,range(len(frgtot.zCMB)),bins=zbins,
        #                                   statistic=lambda values: weighted_avg(values,frgtot,var)).statistic
        #    delmug10simerr = binned_statistic(frgtot.zCMB,range(len(frgtot.zCMB)),bins=zbins,
        #                                      statistic=lambda values: weighted_std(values,frgtot,var)).statistic

        ax.errorbar((zbins[1:]+zbins[:-1])/2.,delmusim,yerr=delmusimerr,fmt='o-',color='C0',label='NIR')
        ax.errorbar((zbins[1:]+zbins[:-1])/2.,delmuoptsim,yerr=delmuoptsimerr,fmt='o-',color='C1',ls='--',label='Opt.+NIR')
        if var == 'DLMAG':
            ax.axhline(0,color='k',lw=2)
            ax.errorbar((zbins[1:]+zbins[:-1])/2.,delmulowz,yerr=delmulowzerr,fmt='^-',color='0.3',ls='-.',label='Pantheon Low-$z$')
            ax.errorbar((zbins[1:]+zbins[:-1])/2.,delmups1,yerr=delmups1err,fmt='^-',color='0.6',ls='-.',label='Pantheon PS1')
            ax.set_ylim([-0.2,0.1])
            #import pdb; pdb.set_trace()
        ax1.legend()

    for ax in [ax1,ax2,ax3]:
        ax.set_xlabel('$z_{CMB}$',fontsize=15)
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    ax1.set_ylabel(r'$\mu - \mu_{sim}$',fontsize=15)
    ax2.set_ylabel(r'$A_V- A_{V,sim}$',fontsize=15)
    ax3.set_ylabel(r'$s_{BV}- s_{BV,sim}$',fontsize=15,labelpad=0)
    plt.savefig('biascor.png')#,bbox_inches='tight',dpi=200)
    
    import pdb; pdb.set_trace()
        
if __name__ == "__main__":
    #main()
    biascor()
