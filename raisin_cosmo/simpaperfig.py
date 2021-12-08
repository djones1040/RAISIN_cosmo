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
import cosmo

# CSP
# median opt stretch:
#    data: 1.02362 +/- 0.127/sqrt(43)
#    sim: 1.06592 +/- 0.0946
# median NIR stretch: 
#    data: 0.9813 +/- 0.115/sqrt(43)
#    sim: 1.11538 +/- 0.113
#
# PS1
# median opt stretch:
#    data: 1.14308 +/- 0.100
#    sim: 0.99538 +/- 0.124
# median NIR stretch: 
#    data: 1.10154 +/- 0.175
#    sim: 1.05538 +/- 0.117
#
# DES
# median opt stretch:
#    data: 1.07111 +/- 0.0929
#    sim: 0.98615 +/- 0.1215
# median NIR stretch: 
#    data: 1.00462 +/- 0.134
#    sim: 1.131 +/- 0.106
#
# weird PS1/DES result.  I think my sims might not be balanced evenly between the two samples?
# and low-z result here is pretty crazy

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

_nirshapesimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN_NIR_SIM_SHAPE',
                          '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN_NIR_SIM_SHAPE',
                          '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN_NIR_SIM_SHAPE']
_nirshapesimfitreslist = [os.path.expandvars(filepath) for filepath in _nirshapesimfitreslist]

#_nirsimfitreslist[1] = '/Users/David/Dropbox/PS1_RAISIN_NIR_SIM_NOSCAT.FITRES.TEXT'
#_nirsimfitreslist[2] = '/Users/David/Dropbox/DES_RAISIN_SIM_NOSCAT.FITRES.TEXT'

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
### ACTUALLY NIR DATA FOR A TEST
#_opticalnirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN_NIR_SHAPECOLOR.FITRES.TEXT',
#                             '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN_NIR_SHAPECOLOR.FITRES.TEXT',
#                             '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN_NIR_SHAPECOLOR.FITRES.TEXT']
_nirshapedatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN_NIR_SHAPE.FITRES.TEXT',
                           '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN_NIR_SHAPE.FITRES.TEXT',
                           '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN_NIR_SHAPE.FITRES.TEXT']
_nirshapedatafitreslist = [os.path.expandvars(filepath) for filepath in _nirshapedatafitreslist]

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
        if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > 0.75 and fropt.STRETCH[fropt.CID == i] < 1.185 and \
           fropt.STRETCHERR[fropt.CID == i] < 0.3: #175:
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


def main(donirshape=True):

    plt.subplots_adjust(wspace=0.0,hspace=0.0)
    
    for nirdatadir,opticalnirdatafitres,opticaldatafitres,nirsimfitres,opticalnirsimfitres,opticalsimfitres,\
        nirshapesimfitres,nirshapedatafitres,nirdatafitres,name,dataversion,idx in zip(
            _outdirs,_opticalnirdatafitreslist,_opticaldatafitreslist,
            _nirsimfitreslist,_opticalnirsimfitreslist,_opticalsimfitreslist,
            _nirshapesimfitreslist,_nirshapedatafitreslist,_nirdatafitreslist,
            ['CSP','PS1','DES'],['CSPDR3','PS1','DES'],range(3)):

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
        try:
            print(glob.glob(os.path.expandvars(f'{opticalnirsimfitres}/{simname}/FITOPT000.FITRES.gz'))[0])
            froptsim = txtobj(glob.glob(os.path.expandvars(f'{opticalnirsimfitres}/{simname}/FITOPT000.FITRES.gz'))[0],
                              fitresheader=True)
            froptsimav = txtobj(glob.glob(os.path.expandvars(f'{opticalnirsimfitres}/{simname}_AVSYS/FITOPT000.FITRES.gz'))[0],
                                fitresheader=True)
            froptsimst = txtobj(glob.glob(os.path.expandvars(f'{opticalnirsimfitres}/{simname}_STSYS/FITOPT000.FITRES.gz'))[0],
                                fitresheader=True)
        except: import pdb; pdb.set_trace()
        iGood = np.array([],dtype=bool)
        for j,i in enumerate(fropt.CID):
            if i in goodcids: iGood = np.append(iGood,True)
            else: iGood = np.append(iGood,False)
        for k in fropt.__dict__.keys():
            fropt.__dict__[k] = fropt.__dict__[k][iGood]
            
        froptsim = apply_all_cuts(froptsim,froptsim,restrict_to_good_list=False)
        froptsimav = apply_all_cuts(froptsimav,froptsimav,restrict_to_good_list=False)
        froptsimst = apply_all_cuts(froptsimst,froptsimst,restrict_to_good_list=False)

        # NIR + shape
        frnirshape = txtobj(nirshapedatafitres,fitresheader=True)
        frnirsimshape = txtobj(glob.glob(os.path.expandvars(f'{nirshapesimfitres}/{simname}/FITOPT000.FITRES.gz'))[0],
                               fitresheader=True)
        frnirsimavshape = txtobj(glob.glob(os.path.expandvars(f'{nirshapesimfitres}/{simname}_AVSYS/FITOPT000.FITRES.gz'))[0],
                                 fitresheader=True)
        frnirsimstshape = txtobj(glob.glob(os.path.expandvars(f'{nirshapesimfitres}/{simname}_STSYS/FITOPT000.FITRES.gz'))[0],
                                 fitresheader=True)

        frnirsimshape = apply_all_cuts(frnirsimshape,froptsim,restrict_to_good_list=False)
        frnirsimavshape = apply_all_cuts(frnirsimavshape,froptsimav,restrict_to_good_list=False)
        frnirsimstshape = apply_all_cuts(frnirsimstshape,froptsimst,restrict_to_good_list=False)
        iGoodShape = (frnirsimshape.STRETCH < 1.165) | (frnirsimshape.STRETCH > 1.185)
        for k in frnirsimshape.__dict__.keys():
            frnirsimshape.__dict__[k] = frnirsimshape.__dict__[k][iGoodShape]
        iGoodShape = (frnirsimavshape.STRETCH < 1.165) | (frnirsimavshape.STRETCH > 1.185)
        for k in frnirsimavshape.__dict__.keys():
            frnirsimavshape.__dict__[k] = frnirsimavshape.__dict__[k][iGoodShape]
        iGoodShape = (frnirsimstshape.STRETCH < 1.165) | (frnirsimstshape.STRETCH > 1.185)
        for k in frnirsimstshape.__dict__.keys():
            frnirsimstshape.__dict__[k] = frnirsimstshape.__dict__[k][iGoodShape]
        
        iGood = np.array([],dtype=bool)
        for j,i in enumerate(frnirshape.CID):
            if i in goodcids: iGood = np.append(iGood,True)
            else: iGood = np.append(iGood,False)
        for k in frnirshape.__dict__.keys():
            frnirshape.__dict__[k] = frnirshape.__dict__[k][iGood]
        #import pdb; pdb.set_trace()

        # NIR no shape
        frnir = txtobj(os.path.expandvars(f"{nirdatadir}/{dataversion}_RAISIN/FITOPT000.FITRES.gz"),fitresheader=True)
        frnirsim = txtobj(glob.glob(os.path.expandvars(f'{nirsimfitres}/{simname}/FITOPT000.FITRES.gz'))[0],
                               fitresheader=True)
        frnirsimav = txtobj(glob.glob(os.path.expandvars(f'{nirsimfitres}/{simname}_AVSYS/FITOPT000.FITRES.gz'))[0],
                            fitresheader=True)
        frnirsimst = txtobj(glob.glob(os.path.expandvars(f'{nirsimfitres}/{simname}_STSYS/FITOPT000.FITRES.gz'))[0],
                            fitresheader=True)

        frnirsim = apply_all_cuts(frnirsim,froptsim,restrict_to_good_list=False)
        frnirsimav = apply_all_cuts(frnirsimav,froptsimav,restrict_to_good_list=False)
        frnirsimst = apply_all_cuts(frnirsimst,froptsimst,restrict_to_good_list=False)
        
        iGood = np.array([],dtype=bool)
        for j,i in enumerate(frnir.CID):
            if i in goodcids: iGood = np.append(iGood,True)
            else: iGood = np.append(iGood,False)
        for k in frnir.__dict__.keys():
            frnir.__dict__[k] = frnir.__dict__[k][iGood]

        # NIR
        hist = ovdatamc.ovhist()
        parser = hist.add_options(usage='')
        options,  args = parser.parse_args()
        hist.options = options
        hist.options.journal = True
        if name != 'CSP': hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.75,1.25),('AV',-0.5,0.45),('SNRMAX1',0,60)]
        else: hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.75,1.25),('AV',-0.5,0.45),('SNRMAX1',0,150)]
        hist.options.nbins = 10
        hist.options.clobber = True
        hist.options.outfile = 'figs/sim_%s.png'%name
        hist.options.histvar = ['MURES'] #,'STRETCH']#,'SNRMAX1']
        hist.main(data=frnir,sim=frnirsim,axes=[ax1],letters=False,ncol=4)
        hist.main(data=frnir,sim=frnirsimav,axes=[ax1],letters=False,simcolor='C0',
                  simlabel='$A_V$ 1$\sigma$',datalabel=None,simls='--',showdata=False,showmeanstd=False,ncol=4)
        hist.main(data=frnir,sim=frnirsimst,axes=[ax1],letters=False,simcolor='C1',
                  simlabel='$s_{BV}$ 1$\sigma$',datalabel=None,simls='-.',showdata=False,showmeanstd=False,ncol=4)

        
        # OPT+NIR
        hist = ovdatamc.ovhist()
        parser = hist.add_options(usage='')
        options,  args = parser.parse_args()
        hist.options = options
        hist.options.journal = True
        if name != 'CSP': hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.75,1.25),('AV',-0.5,0.45),('SNRMAX1',0,60)]
        else: hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.75,1.25),('AV',-0.5,0.45),('SNRMAX1',0,150)]
        hist.options.nbins = 10
        hist.options.clobber = True
        hist.options.outfile = 'figs/sim_%s.png'%name
        if not donirshape:
            hist.options.histvar = ['MURES','AV','STRETCH']#,'SNRMAX1']
            hist.main(data=fropt,sim=froptsim,axes=[ax1,ax2,ax3],letters=False,ncol=4,showlegend=False)
            hist.main(data=fropt,sim=froptsimav,axes=[ax1,ax2,ax3],letters=False,simcolor='C0',
                      simlabel='$A_V$ 1$\sigma$',datalabel=None,simls='--',showdata=False,
                      showmeanstd=False,ncol=4,showlegend=False)
            hist.main(data=fropt,sim=froptsimst,axes=[ax1,ax2,ax3],letters=False,simcolor='C1',
                      simlabel='$s_{BV}$ 1$\sigma$',datalabel=None,simls='-.',showdata=False,
                      showmeanstd=False,ncol=4,showlegend=False)
        else:
            hist.options.histvar = ['AV'] #,'STRETCH']#,'SNRMAX1']
            hist.main(data=fropt,sim=froptsim,axes=[ax2],letters=False,ncol=4,showlegend=False,datalabel=None,ylabel=False)
            hist.main(data=fropt,sim=froptsimav,axes=[ax2],letters=False,simcolor='C0',
                      simlabel='$A_V$ 1$\sigma$',datalabel=None,simls='--',showdata=False,showmeanstd=False,
                      ncol=4,showlegend=False,ylabel=False)
            hist.main(data=fropt,sim=froptsimst,axes=[ax2],letters=False,simcolor='C1',
                      simlabel='$s_{BV}$ 1$\sigma$',datalabel=None,simls='-.',showdata=False,showmeanstd=False,
                      ncol=4,showlegend=False,ylabel=False)


        # NIR + shape
        if donirshape:
            hist = ovdatamc.ovhist()
            parser = hist.add_options(usage='')
            options,  args = parser.parse_args()
            hist.options = options
            hist.options.journal = True
            if name != 'CSP': hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.75,1.25),('AV',-0.5,0.45),('SNRMAX1',0,60)]
            else: hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.75,1.25),('AV',-0.5,0.45),('SNRMAX1',0,150)]
            hist.options.nbins = 10
            hist.options.clobber = True
            hist.options.outfile = 'figs/sim_%s.png'%name
            hist.options.histvar = ['STRETCH']#,'SNRMAX1']
            hist.main(data=frnirshape,sim=frnirsimshape,axes=[ax3],letters=False,ncol=4,showlegend=False,ylabel=False)
            hist.main(data=frnirshape,sim=frnirsimavshape,axes=[ax3],letters=False,simcolor='C0',
                      simlabel='$A_V$ 1$\sigma$',simls='--',showdata=False,showmeanstd=False,ncol=4,datalabel=None,
                      showlegend=False,ylabel=False)
            hist.main(data=frnirshape,sim=frnirsimstshape,axes=[ax3],letters=False,simcolor='C1',
                      simlabel='$s_{BV}$ 1$\sigma$',datalabel=None,simls='-.',showdata=False,showmeanstd=False,ncol=4,
                      showlegend=False,ylabel=False)
            
            
        #import pdb; pdb.set_trace()
        # OPT only
        fropt = txtobj(opticaldatafitres,fitresheader=True)

        print(glob.glob(os.path.expandvars(f'{opticalsimfitres}/{simname}/FITOPT000.FITRES.gz'))[0])
        froptsim = txtobj(glob.glob(os.path.expandvars(f'{opticalsimfitres}/{simname}/FITOPT000.FITRES.gz'))[0],
                          fitresheader=True)

        froptsimav = txtobj(glob.glob(os.path.expandvars(f'{opticalsimfitres}/{simname}_AVSYS/FITOPT000.FITRES.gz'))[0],
                            fitresheader=True)
        froptsimst = txtobj(glob.glob(os.path.expandvars(f'{opticalsimfitres}/{simname}_STSYS/FITOPT000.FITRES.gz'))[0],
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
        if name != 'CSP': hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.75,1.25),('AV',-0.5,0.45),('SNRMAX1',0,60)]
        else: hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.75,1.25),('AV',-0.5,0.45),('SNRMAX1',0,150)]
        #if name != 'CSP': hist.options.cutwin = [('SNRMAX1',0,60)]
        #else: hist.options.cutwin = [('SNRMAX1',0,150)]
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
            if name == 'CSP': ax.set_ylim([0,45])
            else: ax.set_ylim([0,14])
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
            ax1.get_legend().set_bbox_to_anchor([0.949,1.15], transform = ax1.transAxes)
            #ax1.get_legend()._ncol = 4
            #import pdb; pdb.set_trace()
        #if name == 'CSP': import pdb; pdb.set_trace()
        if idx == 0:
            ax1.set_title('NIR')
            ax2.set_title('Opt.$+$NIR')
            ax3.set_title('NIR')
            ax4.set_title('Opt.')
    import pdb; pdb.set_trace()

def biascor(syst=True):
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

        if frfile.endswith('.TEXT'):
            frsim = txtobj(f'{frfile}',fitresheader=True)
        else:
            frsim = txtobj(f'{frfile}/{simname}/FITOPT000.FITRES.gz',fitresheader=True)
            frsimav = txtobj(f'{frfile}/{simname}_AVSYS/FITOPT000.FITRES.gz',fitresheader=True)
            frsimav2 = txtobj(f'{frfile}/{simname}_AVSYS_2/FITOPT000.FITRES.gz',fitresheader=True)
            frsimst = txtobj(f'{frfile}/{simname}_STSYS/FITOPT000.FITRES.gz',fitresheader=True)
        
        froptsim = txtobj(f'{froptfile}/{simname}/FITOPT000.FITRES.gz',fitresheader=True)
        froptsimav = txtobj(f'{froptfile}/{simname}_AVSYS/FITOPT000.FITRES.gz',
                            fitresheader=True)
        froptsimav2 = txtobj(f'{froptfile}/{simname}_AVSYS_2/FITOPT000.FITRES.gz',
                             fitresheader=True)
        froptsimst = txtobj(f'{froptfile}/{simname}_STSYS/FITOPT000.FITRES.gz',
                            fitresheader=True)
        
        frsim = apply_all_cuts(frsim,froptsim,restrict_to_good_list=False)
        frsimav = apply_all_cuts(frsimav,froptsimav,restrict_to_good_list=False)
        frsimav2 = apply_all_cuts(frsimav2,froptsimav2,restrict_to_good_list=False)
        frsimst = apply_all_cuts(frsimst,froptsimst,restrict_to_good_list=False)
        froptsim = apply_all_cuts(froptsim,froptsim,restrict_to_good_list=False)
        froptsimav = apply_all_cuts(froptsimav,froptsimav,restrict_to_good_list=False)
        froptsimav2 = apply_all_cuts(froptsimav2,froptsimav2,restrict_to_good_list=False)
        froptsimst = apply_all_cuts(froptsimst,froptsimst,restrict_to_good_list=False)
        
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
            # let's look at the systematic error variants
            frsimtotav = copy.deepcopy(frsimav)
            froptsimtotav = copy.deepcopy(froptsimav)
            frsimtotav2 = copy.deepcopy(frsimav2)
            froptsimtotav2 = copy.deepcopy(froptsimav2)
            frsimtotst = copy.deepcopy(frsimst)
            froptsimtotst = copy.deepcopy(froptsimst)
        else:
            if name == 'PS1':
                if _nirsimfitreslist[-1].endswith('.TEXT'):
                    froptsimdes = txtobj(f'{_nirsimfitreslist[-1]}',fitresheader=True)
                    froptsimdes = apply_all_cuts(froptsimdes,froptsim,restrict_to_good_list=False)
                    frsimdes = txtobj(f'{_nirsimfitreslist[-1]}',fitresheader=True)
                    frsimdes = apply_all_cuts(frsimdes,froptsimdes,restrict_to_good_list=False)
                else:
                    froptsimdes = txtobj(f'{_nirsimfitreslist[-1]}/DES_RAISIN_SIM/FITOPT000.FITRES.gz',fitresheader=True)
                    froptsimdes = apply_all_cuts(froptsimdes,froptsim,restrict_to_good_list=False)
                    frsimdes = txtobj(f'{_nirsimfitreslist[-1]}/DES_RAISIN_SIM/FITOPT000.FITRES.gz',fitresheader=True)
                    frsimdes = apply_all_cuts(frsimdes,froptsimdes,restrict_to_good_list=False)
                    # syst error variants
                    froptsimdesav = txtobj(f'{_nirsimfitreslist[-1]}/DES_RAISIN_SIM_AVSYS/FITOPT000.FITRES.gz',fitresheader=True)
                    froptsimdesav = apply_all_cuts(froptsimdesav,froptsimav,restrict_to_good_list=False)
                    frsimdesav = txtobj(f'{_nirsimfitreslist[-1]}/DES_RAISIN_SIM_AVSYS/FITOPT000.FITRES.gz',fitresheader=True)
                    frsimdesav = apply_all_cuts(frsimdesav,froptsimdesav,restrict_to_good_list=False)

                    froptsimdesav2 = txtobj(f'{_nirsimfitreslist[-1]}/DES_RAISIN_SIM_AVSYS_2/FITOPT000.FITRES.gz',fitresheader=True)
                    froptsimdesav2 = apply_all_cuts(froptsimdesav2,froptsimav2,restrict_to_good_list=False)
                    frsimdesav2 = txtobj(f'{_nirsimfitreslist[-1]}/DES_RAISIN_SIM_AVSYS_2/FITOPT000.FITRES.gz',fitresheader=True)
                    frsimdesav2 = apply_all_cuts(frsimdesav2,froptsimdesav2,restrict_to_good_list=False)

                    froptsimdesst = txtobj(f'{_nirsimfitreslist[-1]}/DES_RAISIN_SIM_STSYS/FITOPT000.FITRES.gz',fitresheader=True)
                    froptsimdesst = apply_all_cuts(froptsimdesst,froptsimst,restrict_to_good_list=False)
                    frsimdesst = txtobj(f'{_nirsimfitreslist[-1]}/DES_RAISIN_SIM_STSYS/FITOPT000.FITRES.gz',fitresheader=True)
                    frsimdesst = apply_all_cuts(frsimdesst,froptsimdesst,restrict_to_good_list=False)
                
                frgdes = txtobj(f'{_g10fitreslist[-1]}',fitresheader=True)
                frgdes = apply_all_cuts(frgdes,frgdes,restrict_to_good_list=False)
                iFP = frgdes.FITPROB > 1e-3
                for k in frgdes.__dict__.keys():
                    frgdes.__dict__[k] = frgdes.__dict__[k][iFP]
                iZ = frgdes.zHD > 0.1
                for k in frgdes.__dict__.keys():
                    frgdes.__dict__[k] = frgdes.__dict__[k][iZ]
                
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
                
                # syst err
                iRand = sample(range(len(frsimav.CID)),len(frsimdesav.CID))
                iRandOpt = sample(range(len(froptsimav.CID)),len(froptsimdesav.CID))
                for k in frsimtotav.__dict__.keys():
                    frsimtotav.__dict__[k] = np.append(frsimtotav.__dict__[k],frsimav.__dict__[k][iRand])
                for k in froptsimtotav.__dict__.keys():
                    froptsimtotav.__dict__[k] = np.append(froptsimtotav.__dict__[k],froptsimav.__dict__[k][iRandOpt])
                iRand = sample(range(len(frsimav2.CID)),len(frsimdesav2.CID))
                iRandOpt = sample(range(len(froptsimav2.CID)),len(froptsimdesav2.CID))
                for k in frsimtotav2.__dict__.keys():
                    frsimtotav2.__dict__[k] = np.append(frsimtotav2.__dict__[k],frsimav2.__dict__[k][iRand])
                for k in froptsimtotav2.__dict__.keys():
                    froptsimtotav2.__dict__[k] = np.append(froptsimtotav2.__dict__[k],froptsimav2.__dict__[k][iRandOpt])
                iRand = sample(range(len(frsimst.CID)),len(frsimdesst.CID))
                iRandOpt = sample(range(len(froptsimst.CID)),len(froptsimdesst.CID))
                for k in frsimtotst.__dict__.keys():
                    frsimtotst.__dict__[k] = np.append(frsimtotst.__dict__[k],frsimst.__dict__[k][iRand])
                for k in froptsimtotst.__dict__.keys():
                    froptsimtotst.__dict__[k] = np.append(froptsimtotst.__dict__[k],froptsimst.__dict__[k][iRandOpt])
                
            else:

                for k in frgtot.__dict__.keys():
                    frgtot.__dict__[k] = np.append(frgtot.__dict__[k],frg.__dict__[k])
                for k in frsimtot.__dict__.keys():
                    frsimtot.__dict__[k] = np.append(frsimtot.__dict__[k],frsim.__dict__[k])
                for k in froptsimtot.__dict__.keys():
                    froptsimtot.__dict__[k] = np.append(froptsimtot.__dict__[k],froptsim.__dict__[k])
                # syst err
                for k in frsimtotav.__dict__.keys():
                    frsimtotav.__dict__[k] = np.append(frsimtotav.__dict__[k],frsimav.__dict__[k])
                for k in froptsimtotav.__dict__.keys():
                    froptsimtotav.__dict__[k] = np.append(froptsimtotav.__dict__[k],froptsimav.__dict__[k])

                for k in frsimtotav2.__dict__.keys():
                    frsimtotav2.__dict__[k] = np.append(frsimtotav2.__dict__[k],frsimav2.__dict__[k])
                for k in froptsimtotav2.__dict__.keys():
                    froptsimtotav2.__dict__[k] = np.append(froptsimtotav2.__dict__[k],froptsimav2.__dict__[k])

                for k in frsimtotst.__dict__.keys():
                    frsimtotst.__dict__[k] = np.append(frsimtotst.__dict__[k],frsimst.__dict__[k])
                for k in froptsimtotst.__dict__.keys():
                    froptsimtotst.__dict__[k] = np.append(froptsimtotst.__dict__[k],froptsimst.__dict__[k])

                    
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


    zbins = np.linspace(0.01,0.8,15)
    frsimtot.AVERR[:] = 0.01
    frsimtot.STRETCHERR[:] = 0.01

    def get_sigint(resid,residerr):
        chi2_redlist = []
        sigint_testlist = np.arange(0,0.3,0.005)
        for sigint_test in sigint_testlist:
            chi2_redlist += [np.sum(resid**2./(residerr**2.+sigint_test**2.))/float(len(resid)-1)]
        sigint = sigint_testlist[np.abs(np.array(chi2_redlist)-1) == np.min(np.abs(np.array(chi2_redlist)-1))][0]
        return sigint

    # add intrinsic dispersion to reduce influence of outliers
    sigint = get_sigint(frsimtot.DLMAG-cosmo.mu(frsimtot.zHD),frsimtot.DLMAGERR)
    optsigint = get_sigint(froptsimtot.DLMAG-cosmo.mu(froptsimtot.zHD),froptsimtot.DLMAGERR)

    frsimtot.DLMAGERR = np.sqrt(frsimtot.DLMAGERR**2. + sigint**2.)
    froptsimtot.DLMAGERR = np.sqrt(froptsimtot.DLMAGERR**2. + optsigint**2.)
    # syst
    frsimtotav.DLMAGERR = np.sqrt(frsimtotav.DLMAGERR**2. + sigint**2.)
    froptsimtotav.DLMAGERR = np.sqrt(froptsimtotav.DLMAGERR**2. + optsigint**2.)
    frsimtotav2.DLMAGERR = np.sqrt(frsimtotav2.DLMAGERR**2. + sigint**2.)
    froptsimtotav2.DLMAGERR = np.sqrt(froptsimtotav2.DLMAGERR**2. + optsigint**2.)
    frsimtotst.DLMAGERR = np.sqrt(frsimtotst.DLMAGERR**2. + sigint**2.)
    froptsimtotst.DLMAGERR = np.sqrt(froptsimtotst.DLMAGERR**2. + optsigint**2.)

    
    for var,ax in zip(['DLMAG','AV','STRETCH'],[ax1,ax2,ax3]):
        delmusim = binned_statistic(frsimtot.zCMB,range(len(frsimtot.zCMB)),bins=zbins,
                                    statistic=lambda values: weighted_avg(values,frsimtot,var)).statistic
        delmusimerr = binned_statistic(frsimtot.zCMB,range(len(frsimtot.zCMB)),bins=zbins,
                                       statistic=lambda values: weighted_std(values,frsimtot,var)).statistic
        delmusimcount = binned_statistic(frsimtot.zCMB,np.ones(len(frsimtot.DLMAG)),bins=zbins,
                                         statistic='sum').statistic

        delmuoptsim = binned_statistic(froptsimtot.zCMB,range(len(froptsimtot.DLMAG)),bins=zbins,
                                       statistic=lambda values: weighted_avg(values,froptsimtot,var)).statistic
        delmuoptsimerr = binned_statistic(froptsimtot.zCMB,range(len(froptsimtot.DLMAG)),bins=zbins,
                                          statistic=lambda values: weighted_std(values,froptsimtot,var)).statistic
        delmuoptsimcount = binned_statistic(froptsimtot.zCMB,np.ones(len(froptsimtot.DLMAG)),bins=zbins,
                                            statistic='sum').statistic

        
        if var == 'DLMAG':
            delmulowz = binned_statistic(frpanlowz.zCMB,range(len(frpanlowz.zCMB)),bins=zbins,
                                         statistic=lambda values: weighted_avg(values,frpanlowz,var)).statistic
            delmulowzerr = binned_statistic(frpanlowz.zCMB,range(len(frpanlowz.zCMB)),bins=zbins,
                                            statistic=lambda values: weighted_std(values,frpanlowz,var)).statistic
            delmups1 = binned_statistic(frpanps1.zCMB,range(len(frpanps1.zCMB)),bins=zbins,
                                        statistic=lambda values: weighted_avg(values,frpanps1,var)).statistic
            delmups1err = binned_statistic(frpanps1.zCMB,range(len(frpanps1.zCMB)),bins=zbins,
                                           statistic=lambda values: weighted_std(values,frpanps1,var)).statistic

        ax.errorbar((zbins[1:][delmusimcount > 5]+zbins[:-1][delmusimcount > 5])/2.,
                    delmusim[delmusimcount > 5],yerr=delmusimerr[delmusimcount > 5],
                    fmt='o-',color='C0',label='NIR')
        ax.errorbar((zbins[1:][delmuoptsimcount > 5]+zbins[:-1][delmuoptsimcount > 5])/2.,
                    delmuoptsim[delmuoptsimcount > 5],yerr=delmuoptsimerr[delmuoptsimcount > 5],
                    fmt='o-',color='C1',ls='--',label='Opt.+NIR')
        # syst
        if syst:
            # syst
            delmusimav = binned_statistic(frsimtotav.zCMB,range(len(frsimtotav.DLMAG)),bins=zbins,
                                          statistic=lambda values: weighted_avg(values,frsimtotav,var)).statistic
            delmusimaverr = binned_statistic(frsimtotav.zCMB,range(len(frsimtotav.DLMAG)),bins=zbins,
                                             statistic=lambda values: weighted_std(values,frsimtotav,var)).statistic
            delmusimavcount = binned_statistic(frsimtotav.zCMB,np.ones(len(frsimtotav.DLMAG)),bins=zbins,
                                               statistic='sum').statistic
        
            delmuoptsimav = binned_statistic(froptsimtotav.zCMB,range(len(froptsimtotav.DLMAG)),bins=zbins,
                                             statistic=lambda values: weighted_avg(values,froptsimtotav,var)).statistic
            delmuoptsimaverr = binned_statistic(froptsimtotav.zCMB,range(len(froptsimtotav.DLMAG)),bins=zbins,
                                                statistic=lambda values: weighted_std(values,froptsimtotav,var)).statistic
            delmuoptsimavcount = binned_statistic(froptsimtotav.zCMB,np.ones(len(froptsimtotav.DLMAG)),bins=zbins,
                                                  statistic='sum').statistic

            delmusimav2 = binned_statistic(frsimtotav2.zCMB,range(len(frsimtotav2.DLMAG)),bins=zbins,
                                          statistic=lambda values: weighted_avg(values,frsimtotav2,var)).statistic
            delmusimav2err = binned_statistic(frsimtotav2.zCMB,range(len(frsimtotav2.DLMAG)),bins=zbins,
                                             statistic=lambda values: weighted_std(values,frsimtotav2,var)).statistic
            delmusimav2count = binned_statistic(frsimtotav2.zCMB,np.ones(len(frsimtotav2.DLMAG)),bins=zbins,
                                               statistic='sum').statistic
        
            delmuoptsimav2 = binned_statistic(froptsimtotav2.zCMB,range(len(froptsimtotav2.DLMAG)),bins=zbins,
                                             statistic=lambda values: weighted_avg(values,froptsimtotav2,var)).statistic
            delmuoptsimav2err = binned_statistic(froptsimtotav2.zCMB,range(len(froptsimtotav2.DLMAG)),bins=zbins,
                                                statistic=lambda values: weighted_std(values,froptsimtotav2,var)).statistic
            delmuoptsimav2count = binned_statistic(froptsimtotav2.zCMB,np.ones(len(froptsimtotav2.DLMAG)),bins=zbins,
                                                  statistic='sum').statistic

            delmusimst = binned_statistic(frsimtotav.zCMB,range(len(frsimtotav.DLMAG)),bins=zbins,
                                          statistic=lambda values: weighted_avg(values,frsimtotav,var)).statistic
            delmusimsterr = binned_statistic(frsimtotav.zCMB,range(len(frsimtotav.DLMAG)),bins=zbins,
                                             statistic=lambda values: weighted_std(values,frsimtotav,var)).statistic
            delmusimstcount = binned_statistic(frsimtotav.zCMB,np.ones(len(frsimtotav.DLMAG)),bins=zbins,
                                               statistic='sum').statistic
        
            delmuoptsimst = binned_statistic(froptsimtotst.zCMB,range(len(froptsimtotst.DLMAG)),bins=zbins,
                                             statistic=lambda values: weighted_avg(values,froptsimtotst,var)).statistic
            delmuoptsimsterr = binned_statistic(froptsimtotst.zCMB,range(len(froptsimtotst.DLMAG)),bins=zbins,
                                                statistic=lambda values: weighted_std(values,froptsimtotst,var)).statistic
            delmuoptsimstcount = binned_statistic(froptsimtotst.zCMB,np.ones(len(froptsimtotst.DLMAG)),bins=zbins,
                                                  statistic='sum').statistic

            
            sys_val = None #'av2'
            if sys_val == 'av':
                ax.errorbar((zbins[1:][delmusimavcount > 5]+zbins[:-1][delmusimavcount > 5])/2.,
                            delmusimav[delmusimavcount > 5],yerr=delmusimaverr[delmusimavcount > 5],
                            fmt='o-',color='C2',ls='--',label='AV Sys (NIR)')
                ax.errorbar((zbins[1:][delmuoptsimavcount > 5]+zbins[:-1][delmuoptsimavcount > 5])/2.,
                            delmuoptsimav[delmuoptsimavcount > 5],yerr=delmuoptsimaverr[delmuoptsimavcount > 5],
                            fmt='o-',color='C3',ls='--',label='AV Sys (Opt+NIR)')
            elif sys_val == 'av2':
                ax.errorbar((zbins[1:][delmusimav2count > 5]+zbins[:-1][delmusimav2count > 5])/2.,
                            delmusimav2[delmusimav2count > 5],yerr=delmusimav2err[delmusimav2count > 5],
                            fmt='o-',color='C2',ls='--',label='AV2 Sys (NIR)')
                ax.errorbar((zbins[1:][delmuoptsimav2count > 5]+zbins[:-1][delmuoptsimav2count > 5])/2.,
                            delmuoptsimav2[delmuoptsimav2count > 5],yerr=delmuoptsimav2err[delmuoptsimav2count > 5],
                            fmt='o-',color='C3',ls='--',label='AV2 Sys (Opt+NIR)')
            elif sys_val == 'st':
                ax.errorbar((zbins[1:][delmusimstcount > 5]+zbins[:-1][delmusimstcount > 5])/2.,
                            delmusimst[delmusimstcount > 5],yerr=delmusimsterr[delmusimstcount > 5],
                            fmt='o-',color='C2',ls='--',label='ST Sys (NIR)')
                ax.errorbar((zbins[1:][delmuoptsimstcount > 5]+zbins[:-1][delmuoptsimstcount > 5])/2.,
                            delmuoptsimst[delmuoptsimstcount > 5],yerr=delmuoptsimsterr[delmuoptsimstcount > 5],
                            fmt='o-',color='C3',ls='--',label='ST Sys (Opt+NIR)')

            
        if var == 'DLMAG':
            ax.axhline(0,color='k',lw=2)
            ax.errorbar((zbins[1:]+zbins[:-1])/2.,delmulowz,yerr=delmulowzerr,fmt='^-',color='0.3',ls='-.',label='Pantheon Low-$z$')
            ax.errorbar((zbins[1:]+zbins[:-1])/2.,delmups1,yerr=delmups1err,fmt='^-',color='0.6',ls='-.',label='Pantheon PS1')
            ax.set_ylim([-0.2,0.12])

#        import pdb; pdb.set_trace()
        ax1.legend()

    for ax in [ax1,ax2,ax3]:
        ax.set_xlabel('$z_{CMB}$',fontsize=15)
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    ax1.set_ylabel(r'$\mu_{fit} - \mu_{sim}$ (mag)',fontsize=15)
    ax2.set_ylabel(r'$A_{V,fit}- A_{V,sim}$ (mag)',fontsize=15)
    ax3.set_ylabel(r'$s_{BV,fit}- s_{BV,sim}$',fontsize=15,labelpad=0)
    plt.savefig('biascor.png')#,bbox_inches='tight',dpi=200)
    
    import pdb; pdb.set_trace()
        
if __name__ == "__main__":
    #main()
    biascor()
