#!/usr/bin/env python
# D. Jones - 9/29/21
"""
bias corrections that take the true stretch
but otherwise use the NIR distances
"""

import numpy as np
from astropy.io import fits
import glob
import os
from txtobj import txtobj
from scipy.stats import binned_statistic

_nmltmpl = """
  &SNLCINP
	 PRIVATE_DATA_PATH = '$RAISIN_ROOT/cosmo/data/Photometry'

	 VERSION_PHOTOMETRY = '<data_version>'
	 KCOR_FILE		   = '$RAISIN_ROOT/cosmo/kcor/<kcor>'

	 NFIT_ITERATION = 3
	 INTERP_OPT		= 1

	 SNTABLE_LIST = 'FITRES LCPLOT(text:key)'
	 TEXTFILE_PREFIX  = '<outfile>'
	 
	 LDMP_SNFAIL = T
	 USE_MWCOR = F
     USE_MINOS = F

	 H0_REF	  = 70.0
	 OLAM_REF =	 0.70
	 OMAT_REF =	 0.30
	 W0_REF	  = -1.00

	 SNCID_LIST	   =  0
	 CUTWIN_CID	   =  0, 20000000
	 SNCCID_LIST   =  <cidlist>
	 SNCCID_IGNORE =  

	 cutwin_redshift   = 0.001, 2.0
	 cutwin_Nepoch	  =	 1

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

	 FITMODEL_NAME	= '$RAISIN_ROOT/cosmo/snoopy.B18'
	 OPT_PRIOR_AV = 0
    
	 PRIOR_MJDSIG		 = 5.0
	 PRIOR_LUMIPAR_RANGE = -5.0, 5.0
     INIVAL_GRIDSEARCH_COLOR = -1.0, 1.0, 0.05
     !INIVAL_SHAPE = <inival_st>
     !INISTP_SHAPE = 0.0
     !INIVAL_AV = <inival_av>
     !INISTP_AV = 0.0
	 !INISTP_PEAKMJD = 0.0
     !INIVAL_DLMAG = <inival_dlmag>
     !INISTP_DLMAG = 0.0

	 INIVAL_RV = 1.518

	 !OPT_COVAR = 1
	 OPT_XTMW_ERR = 1
	 OPT_COVAR_FLUX = 0
	 TREST_REJECT  = -15.0, 45.0
	 NGRID_PDF	   = 0

	 FUDGEALL_ITER1_MAXFRAC = 0.02
	 FILTLIST_FIT = '<filtlist>'

  &END
"""


def apply_all_cuts(fr,fropt,restrict_to_good_list=False):

    # AV
    iGoodAV = np.zeros(len(fr.CID),dtype=bool)
    for j,i in enumerate(fr.CID):
        if i in fropt.CID and fropt.AV[fropt.CID == i] < 0.3*fropt.RV[fropt.CID == i]:
            iGoodAV[j] = True

    # reasonable stretch
    iGoodSt = np.zeros(len(fr.CID),dtype=bool)
    for j,i in enumerate(fr.CID):
        if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > 0.75 and fropt.STRETCH[fropt.CID == i] < 1.185: #175:
            iGoodSt[j] = True

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


class biascor:
    def __init__(self):
        self.band = 'Y'
        self.versions = ['CSP_RAISIN_SIM','PS1_RAISIN_SIM','DES_RAISIN_SIM']
        self.CIDS = [np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                     np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                     np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)]
        self.kcors = ['kcor_CSPDR3_BD17.fits','kcor_PS1MD_NIR.fits','kcor_DES_NIR.fits']
        self.filtlist = ['YyJjH','JH','JH']
    
    def run_snana(self,cidlist,version,kcor,filtlist,outfile,
                  inival_st=None,inival_av=None,inival_dlmag=None,clobber=False):

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

        
        for cids,kcor,filtlist,version in zip(
                self.CIDS,self.kcors,self.filtlist,self.versions):

            hdu = fits.open(os.path.expandvars(f'$SNDATA_ROOT/SCRATCH/{version}/{version}_HEAD.FITS'))
            cids = hdu[1].data['SNID']
            stretch = hdu[1].data['SIM_STRETCH']
            with open(f'{version}_biascor_stretch_combined.FITRES.TEXT','w') as fout:
                print('# CID zCMB STRETCH STRETCHERR AV AVERR DLMAG DLMAGERR SIM_DLMAG',file=fout)
                for c,s in zip(cids[0:1000],stretch[0:1000]):
                
                    outfile = self.run_snana([c],version,kcor,filtlist,f"biascor_stretch/{version}_biascor_stretch_{c}",
                                             clobber=False,inival_av=0,inival_st=s)
                
                    fr = txtobj(f"{outfile}.FITRES.TEXT",fitresheader=True)
                    print(f"{fr.CID[0]:.0f} {fr.zHD[0]:.5f} {fr.STRETCH[0]:.5f} {fr.STRETCHERR[0]:.5f} {fr.AV[0]:.5f} {fr.AVERR[0]:.5f} {fr.DLMAG[0]:.5f} {fr.DLMAGERR[0]:.5f} {fr.SIM_DLMAG[0]:.5f}",file=fout)
                    #import pdb; pdb.set_trace()
    
    def biascor(self):
        # plot up the biascors and write to files

        zbins_lowz = np.linspace(0,0.1,7)
        zbins_highz = np.linspace(0.1,0.7,10)

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

        # low-z
        frsimtot = txtobj(f'{self.versions[0]}_biascor_stretch_combined.FITRES.TEXT')
        froptsimtot = txtobj('output/fit_all/CSP_RAISIN_OPTNIR_SIM/CSP_RAISIN_SIM/FITOPT000.FITRES.gz',fitresheader=True)
        frsimtot = apply_all_cuts(frsimtot,froptsimtot)
        frsimtot.DLMAGERR = np.sqrt(frsimtot.DLMAGERR**2.+0.1**2.)
        delmusim = binned_statistic(frsimtot.zCMB,range(len(frsimtot.zCMB)),bins=zbins_lowz,
                                    statistic=lambda values: weighted_avg(values,frsimtot)).statistic
        delmusimerr = binned_statistic(frsimtot.zCMB,range(len(frsimtot.zCMB)),bins=zbins_lowz,
                                       statistic=lambda values: weighted_std(values,frsimtot)).statistic
        with open('CSP_biascor_stretch.txt','w') as fout:
            for z,dm,dme in zip((zbins_lowz[1:]+zbins_lowz[:-1])/2.,delmusim,delmusimerr):
                print(f"{z:.3f} {dm:.3f} {dme:.3f}",file=fout)
        
        # PS1
        frsimtot = txtobj(f'{self.versions[1]}_biascor_stretch_combined.FITRES.TEXT')
        froptsimtot = txtobj('output/fit_all/PS1_RAISIN_OPTNIR_SIM/PS1_RAISIN_SIM/FITOPT000.FITRES.gz',fitresheader=True)
        frsimtot = apply_all_cuts(frsimtot,froptsimtot)
        frsimtot.DLMAGERR = np.sqrt(frsimtot.DLMAGERR**2.+0.1**2.)
        delmusim = binned_statistic(frsimtot.zCMB,range(len(frsimtot.zCMB)),bins=zbins_highz,
                                    statistic=lambda values: weighted_avg(values,frsimtot)).statistic
        delmusimerr = binned_statistic(frsimtot.zCMB,range(len(frsimtot.zCMB)),bins=zbins_highz,
                                       statistic=lambda values: weighted_std(values,frsimtot)).statistic
        with open('PS1_biascor_stretch.txt','w') as fout:
            for z,dm,dme in zip((zbins_highz[1:]+zbins_highz[:-1])/2.,delmusim,delmusimerr):
                print(f"{z:.3f} {dm:.3f} {dme:.3f}",file=fout)

        # DES
        frsimtot = txtobj(f'{self.versions[1]}_biascor_stretch_combined.FITRES.TEXT')
        froptsimtot = txtobj('output/fit_all/DES_RAISIN_OPTNIR_SIM/DES_RAISIN_SIM/FITOPT000.FITRES.gz',fitresheader=True)
        frsimtot = apply_all_cuts(frsimtot,froptsimtot)
        frsimtot.DLMAGERR = np.sqrt(frsimtot.DLMAGERR**2.+0.1**2.)
        delmusim = binned_statistic(frsimtot.zCMB,range(len(frsimtot.zCMB)),bins=zbins_highz,
                                    statistic=lambda values: weighted_avg(values,frsimtot)).statistic
        delmusimerr = binned_statistic(frsimtot.zCMB,range(len(frsimtot.zCMB)),bins=zbins_highz,
                                       statistic=lambda values: weighted_std(values,frsimtot)).statistic
        with open('DES_biascor_stretch.txt','w') as fout:
            for z,dm,dme in zip((zbins_highz[1:]+zbins_highz[:-1])/2.,delmusim,delmusimerr):
                print(f"{z:.3f} {dm:.3f} {dme:.3f}",file=fout)

        
if __name__ == "__main__":
    bc = biascor()
    bc.main()
    #bc.biascor()
