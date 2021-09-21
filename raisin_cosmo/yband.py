#!/usr/bin/env python
# D. Jones - 9/11/21
# take a look at some Y-band fitting options in the low-z sample

from txtobj import txtobj
import pylab as plt
from matplotlib.gridspec import GridSpec
plt.ion()
import numpy as np
import cosmo
import os
from scipy.stats import binned_statistic
import glob
from astropy.stats import sigma_clipped_stats
from scipy.interpolate import interp1d

goodcids = np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)

def Yfits():
    fre = txtobj('output/fit_nir/CSP_RAISIN_onlyY_early.FITRES.TEXT',fitresheader=True)
    frl = txtobj('output/fit_nir/CSP_RAISIN_onlyY_late.FITRES.TEXT',fitresheader=True)

    idx = np.array([],dtype=int)
    for j,i in enumerate(fre.CID):
        if i in goodcids:
            idx = np.append(idx,j)
    for k in fre.__dict__.keys():
        fre.__dict__[k] = fre.__dict__[k][idx]
    idx = np.array([],dtype=int)
    for j,i in enumerate(frl.CID):
        if i in goodcids:
            idx = np.append(idx,j)
    for k in frl.__dict__.keys():
        frl.__dict__[k] = frl.__dict__[k][idx]
    
    
    iFPe = fre.FITPROB > -1 #1e-6
    iFPl = frl.FITPROB > -1 #1e-6
    plt.errorbar(fre.zHD[iFPe],fre.DLMAG[iFPe]-cosmo.mu(fre.zHD[iFPe]),fmt='o')
    plt.errorbar(frl.zHD[iFPl],frl.DLMAG[iFPl]-cosmo.mu(frl.zHD[iFPl]),fmt='o')

    print(np.median(fre.DLMAG[iFPe]-cosmo.mu(fre.zHD[iFPe])))
    print(np.median(frl.DLMAG[iFPl]-cosmo.mu(frl.zHD[iFPl])))
    
    import pdb; pdb.set_trace()

def DESY():

    fr = txtobj('output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT',fitresheader=True)
    DES_late_cids = ['DES15C3odz',
                     'DES16C1cim', 'DES16S1bno', 'DES16C2cva', 'DES15E2nlz',
                     'DES15E2uc', 'DES15X2kvt', 'DES16E2clk', 'DES16E2cqq']
    DES_early_cids = ['DES15X2nkz', 'DES16S1agd', 'DES15E2mhy', 'DES16S2afz',
                      'DES16C3cmy', 'DES16X3cry', 'DES16E1dcx', 'DES16X3zd',
                      'DES15X2mey']
    sbvl = np.array([]); avl = np.array([])
    for i in DES_late_cids:
        sbvl = np.append(sbvl,fr.STRETCH[fr.CID == i])
        avl = np.append(avl,fr.AV[fr.CID == i])
    sbve = np.array([]); ave = np.array([])
    for i in DES_early_cids:
        sbve = np.append(sbve,fr.STRETCH[fr.CID == i])
        ave = np.append(ave,fr.AV[fr.CID == i])
    print(np.median(sbve),np.median(sbvl))
    print(np.median(ave),np.median(avl))
    import pdb; pdb.set_trace()

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
    
class YbandModel:
    def __init__(self):
        self.band = 'Y'
        self.versions = ['CSPDR3_RAISIN','PS1_RAISIN','DES_RAISIN']
        self.CIDS = [np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                     np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                     np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)]
        self.kcors = ['kcor_CSPDR3_BD17.fits','kcor_PS1MD_NIR.fits','kcor_DES_NIR.fits']
        self.filtlist = ['BomngriYyJjH','grizJH','grizJH']

        # 40.80119 -> 40.7813
        
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
        #os.system('rm tmp.nml')

        return outfile
        
    def main(self):

        fig = plt.figure()#constrained_layout=True)
        plt.subplots_adjust(top=0.95)
        gs = GridSpec(3, 3, figure=fig)
        gs.update(wspace=0.0, hspace=0.4)

        #plt.tick_params(left)

        axmain = fig.add_subplot(gs[1:, :])
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])
                
        # 1. get CID list
        # 2. run opt+NIR SNANA fits
        phase,dataflux,datafluxerr,modelflux,survey,cidout = np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
        stretchout,colorout,zout = np.array([]),np.array([]),np.array([])
        for cids,kcor,filtlist,version,ax,title in zip(
                self.CIDS,self.kcors,self.filtlist,self.versions,[ax1,ax2,ax3],['CSP','PS1 (RAISIN1)','DES (RAISIN2)']):
            #if 'CSP' in version: continue
            outfile = self.run_snana(cids,version,kcor,filtlist,f"{self.band.lower()}band/{version}_optnir",clobber=True)
            #outfile = f"{version}_optnir"

            # 3. figure out which obs bands -> Y (or whatever band)
            lcm = txtobj(f"{outfile}.LCPLOT.TEXT",fitresheader=True,rowprfx='OBS')
            frm = txtobj(f"{outfile}.FITRES.TEXT",fitresheader=True)

            cid_filtdict = {}
            for c in np.unique(lcm.CID):
                # move on to the next SN if no Y band
                if self.band not in lcm.BAND_REST[lcm.CID == c]: continue
                if frm.AV[frm.CID == c][0] > 1.0 or frm.STRETCH[frm.CID == c][0] < 0.8 or frm.STRETCH[frm.CID == c][0] > 1.3:
                    continue
                # move on to the next SN if it fails cuts
                
                bands_rest = lcm.BAND_REST[lcm.CID == c]
                iRest = (bands_rest == self.band) #| (bands_rest == 'J')
                bands_to_avoid = np.unique(lcm.BAND[lcm.CID == c][iRest])
                #if len(bands_to_avoid) == 1:
                cid_filtdict[c] = bands_to_avoid #[0]
                #else: raise RuntimeError("I'm confused")
                
            
                # 4. run opt+NIR SNANA fits w/o that band (probably do one SN at a time)
                pardict = {}
                outfile = self.run_snana(
                    [c],version,kcor,"".join([f for f in filtlist if f not in bands_to_avoid]),
                    f"{self.band.lower()}band/{version}_optnir_single_{c}",clobber=True)
                fr = txtobj(f"{outfile}.FITRES.TEXT",fitresheader=True)
                try: pardict[c] = (fr.DLMAG[0],fr.STRETCH[0],fr.AV[0],fr.zHEL[0])
                except: continue
                #print("".join([f for f in filtlist if f not in bands_to_avoid]))
                #import pdb; pdb.set_trace()
                # 5. run SNANA fits w/ Y-band only and no DLMAG (for other code)
                outfile = self.run_snana(
                    [c],version,kcor,"".join(bands_to_avoid), #n+"_modeladjust"
                    f"{self.band.lower()}band/{version}_{self.band}_single_{c}_nodlmag",
                    inival_st=fr.STRETCH[0],inival_av=fr.AV[0],clobber=False) #inival_dlmag=fr.DLMAG[0]
                frmt = txtobj(f"{outfile}.FITRES.TEXT",fitresheader=True)
                if fr.STRETCH[0] <1.05: continue
                if fr.AV[0] > 0.4: continue
                
                # 5. run SNANA fits w/ Y-band only
                outfile = self.run_snana(
                    [c],version,kcor,"".join(bands_to_avoid),
                    f"{self.band.lower()}band/{version}_{self.band}_single_{c}",
                    inival_st=fr.STRETCH[0],inival_av=fr.AV[0],inival_dlmag=fr.DLMAG[0],clobber=True)
                
                # 6. get model from the LCPLOT files for each SN
                lc = txtobj(f"{outfile}.LCPLOT.TEXT",fitresheader=True,rowprfx='OBS')
                frs = txtobj(f"{outfile}.FITRES.TEXT",fitresheader=True)

                for mjd,df,dfe in zip(
                        lc.MJD[lc.DATAFLAG == 1],lc.FLUXCAL[lc.DATAFLAG == 1],lc.FLUXCAL_ERR[lc.DATAFLAG == 1]):
                    phase = np.append(phase,(mjd-frs.PKMJD[0])/(1+fr.zHEL[0]))
                    dataflux = np.append(dataflux,df)
                    datafluxerr = np.append(datafluxerr,dfe)
                    mf = np.interp([mjd-frs.PKMJD[0]],lc.Tobs[lc.DATAFLAG == 0],lc.FLUXCAL[lc.DATAFLAG == 0])[0]
                    # scale the modelflux to match the amplitude of the previous fit
                    modelflux = np.append(modelflux,mf) #*10**(0.4*fr.DLMAG[0]-frs.DLMAG[0]))
                    survey = np.append(survey,version)
                    cidout = np.append(cidout,c)
                    stretchout = np.append(stretchout,fr.STRETCH[0])
                    colorout = np.append(colorout,fr.AV[0])
                    zout = np.append(zout,fr.zHD[0])
                #if c == 'DES15E2nlz': import pdb; pdb.set_trace()
            
            datamag = -2.5*np.log10(dataflux)+27.5
            datamagerr = 1.086*datafluxerr/dataflux
            modelmag = -2.5*np.log10(modelflux)+27.5
            ax.errorbar(
                phase[survey == version],datamag[survey == version]-modelmag[survey == version],
                yerr=datamagerr[survey == version],fmt='.')
            ax.set_title(title)
            ax.axhline(0,color='k')
            ax.set_xlabel('Phase')
            #ax.legend()
            ax.set_ylim([-0.5,0.5])
            #import pdb; pdb.set_trace()
            
            N = 1 #int(len(datamag[survey == version])/15) # arbitrary
            phase_sort = np.sort(phase[survey == version])
            residmag = datamag[survey == version]-modelmag[survey == version]
            residmag_sort = residmag[np.argsort(phase[survey == version])]

            rolling_mean = np.convolve(residmag_sort, np.ones((N,))/N, mode='same')
            #axmain.plot(phase_sort,rolling_mean,zorder=100,lw=4)
            with open(f'modeladjust_{version}.txt','w') as fout:
                for p,rm in zip(phase_sort,rolling_mean):
                    print(f"{p:.3f} {rm:.3f}",file=fout)

        # 8. run non-parametric fit to get suggested model adjustments, weighting each SN equally somehow
        #    or a rolling median?
        N = int(len(datamag)/15) # arbitrary
        phase_sort = np.sort(phase)
        residmag = datamag-modelmag
        residmag_sort = residmag[np.argsort(phase)]

        ax1.set_ylabel(f"${self.band}$ data-model")
        #axmain.plot(phase_sort,np.convolve(residmag_sort, np.ones((N,))/N, mode='same'),zorder=100,lw=4)

            #import pdb; pdb.set_trace()
        # 7. compare data to model
        datamag = -2.5*np.log10(dataflux)+27.5
        datamagerr = 1.086*datafluxerr/dataflux
        modelmag = -2.5*np.log10(modelflux)+27.5
        #for version in self.versions:
        #    plt.errorbar(
        #        phase[survey == version],datamag[survey == version]-modelmag[survey == version],
        #        yerr=datamagerr[survey == version],fmt='o',label=version)
        axmain.axhline(0,color='k')
        axmain.set_xlabel('Phase')
        axmain.set_ylabel(f"${self.band}$ data-model")
        axmain.legend()
        
        # 8. run non-parametric fit to get suggested model adjustments, weighting each SN equally somehow
        #    or a rolling median?
        N = int(len(datamag)/15) # arbitrary
        phase_sort = np.sort(phase)
        residmag = datamag-modelmag
        residmag_sort = residmag[np.argsort(phase)]
        survey_sort = survey[np.argsort(phase)]
        
        N = 50
        axmain.plot(phase_sort[survey_sort == 'CSPDR3_RAISIN'],
                    np.convolve(residmag_sort[survey_sort == 'CSPDR3_RAISIN'], np.ones((N,))/N, mode='same'),
                    zorder=100,lw=4,label='CSP')
        N = 10
        axmain.plot(phase_sort[survey_sort != 'CSPDR3_RAISIN'],
                    np.convolve(residmag_sort[survey_sort != 'CSPDR3_RAISIN'],np.ones((N,))/N, mode='same'),
                    zorder=100,lw=4,label='RAISIN')
        ax2.yaxis.set_ticklabels([])
        ax3.yaxis.tick_right()
        with open('modeladjust_highz.txt','w') as fout:
            for p,m in zip(
                    phase_sort[survey_sort != 'CSPDR3_RAISIN'],
                    np.convolve(residmag_sort[survey_sort != 'CSPDR3_RAISIN'],np.ones((N,))/N, mode='same')):
                print(f"{p:.4f} {m:.4f}",file=fout)                
        
        # 9. now get the "final" model adjustment by taking a straight average of the low- and high-z results
        #    though we might want to test using only the high-z results to see if Kaisey's trend disappears
        phase_bins = np.arange(-20,50,2.0)


        phase_sort = np.sort(phase[survey != 'CSPDR3_RAISIN'])
        residmag = datamag[survey != 'CSPDR3_RAISIN']-modelmag[survey != 'CSPDR3_RAISIN']
        residmag_sort = residmag[np.argsort(phase[survey != 'CSPDR3_RAISIN'])]
        def weighted_avg(idx):
            if not len(idx): return np.nan
            average = np.average(residmag_sort[idx], weights=1/datamagerr[idx]**2.)
            return average
        def sigma_clipped_avg(idx):
            if not len(idx): return np.nan
            average = sigma_clipped_stats(residmag_sort[idx])[0]
            return average

        mean_highz = binned_statistic(phase_sort,range(len(residmag_sort)),bins=phase_bins,statistic=sigma_clipped_avg).statistic
        axmain.plot((phase_bins[1:]+phase_bins[:-1])/2.,mean_highz,label='sigma_clipped mean (RAISIN)')
        #import pdb; pdb.set_trace()
        #return
        phase_sort = np.sort(phase[survey == 'CSPDR3_RAISIN'])
        residmag = datamag[survey == 'CSPDR3_RAISIN']-modelmag[survey == 'CSPDR3_RAISIN']
        residmag_sort = residmag[np.argsort(phase[survey == 'CSPDR3_RAISIN'])]
        def weighted_avg(idx):
            if not len(idx): return np.nan
            average = np.average(residmag_sort[idx], weights=1/datamagerr[idx]**2.)
            return average
        def sigma_clipped_avg(idx):
            if not len(idx): return np.nan
            average = sigma_clipped_stats(residmag_sort[idx])[0]
            return average

        mean_lowz = binned_statistic(phase_sort,range(len(residmag_sort)),bins=phase_bins,statistic=sigma_clipped_avg).statistic

        axmain.legend()
        
        with open(f'modeladjust_final.txt','w') as fout:
            for p1,p2,rmh,rml in zip(phase_bins[:-1],phase_bins[1:],mean_highz,mean_lowz):
                print(f"{(p1+p2)/2.:.4f} {np.nanmean([rmh]):.4f}",file=fout)
        plt.savefig('raisin_modeladjust.png',dpi=200)
        
        import pdb; pdb.set_trace()

    def adjust_data(self):
        ###
        # we can't easily adjust the model
        # but we can instead adjust the data to eliminate data-model residuals as a function of phase
        ###

        phase,resid = np.loadtxt('modeladjust_final.txt',unpack=True) #'modeladjust_final.txt',unpack=True)
        for version in ['CSPDR3_RAISIN','PS1_RAISIN','DES_RAISIN'][1:]:
            lcm = txtobj(f"yband/{version}_optnir.LCPLOT.TEXT",fitresheader=True,rowprfx='OBS')
            frm = txtobj(f"yband/{version}_optnir.FITRES.TEXT",fitresheader=True)

            newdir = f"data/Photometry/{version}_modeladjust"
            if not os.path.exists(newdir):
                os.makedirs(newdir)
            files = glob.glob(f"data/Photometry/{version}/*")
            for f in files:
                varnames = None
                with open(f) as fin, open(f.replace(version,f"{version}_modeladjust"),"w") as fout:
                    for line in fin:
                        line = line.replace('\n','')
                        if line.startswith('VARLIST:'):
                            varnames = np.array(line.split())
                            print(line,file=fout)
                        elif line.startswith('SNID:'):
                            snid = line.split()[1]
                            print(line,file=fout)
                        elif line.startswith('PEAKMJD'):
                            #try: pkmjd = frm.PKMJD[frm.CID == snid][0]
                            pkmjd = float(line.split()[1])
                            print(line,file=fout)
                        elif line.startswith('REDSHIFT_HELIO'):
                            zhel = float(line.split()[1])
                            print(line,file=fout)
                        elif line.startswith('OBS:'):
                            lineparts = np.array(line.split())
                            flt = lineparts[(varnames == 'FLT') | (varnames == 'BAND')][0]
                            band_rest = np.unique(lcm.BAND_REST[(lcm.CID == snid) & (lcm.BAND == flt) & (lcm.DATAFLAG == 1)])
                            if not len(band_rest):
                                print(line,file=fout)
                            elif band_rest[0] != 'Y': #flt != 'Y' and flt != 'y':
                                print(line,file=fout)
                            else:
                                ifluxcal = np.where(varnames == 'FLUXCAL')[0]
                                mjd = float(lineparts[varnames == 'MJD'][0])
                                singlephase = (mjd-pkmjd)/(1+zhel)
                                magadjust = np.interp(singlephase,phase[resid == resid],resid[resid == resid])
                                fluxcal_new = float(lineparts[ifluxcal])*10**(0.4*magadjust) # data becomes brighter to match the model
                                if fluxcal_new != fluxcal_new:
                                    import pdb; pdb.set_trace()
                                lineparts[ifluxcal] = fluxcal_new
                                print(" ".join(lineparts),file=fout)
                                #if snid == 'PScA470110': import pdb; pdb.set_trace()
                        else:
                            print(line,file=fout)
                #import pdb; pdb.set_trace()
            os.chdir(f"data/Photometry/{version}_modeladjust")
            os.system(f'ls *dat > {version}_modeladjust.LIST')
            os.system(f'ls *DAT >> {version}_modeladjust.LIST')
            os.chdir("../../../")
            with open(f"data/Photometry/{version}_modeladjust/{version}_modeladjust.README","w") as fout:
                print("""DOCUMENTATION_START:
blah!
DOCUMENTATION_END:""",file=fout)

    def adjust_data_bayesn(self):
        ###
        # we can't easily adjust the model
        # but we can instead adjust the data to eliminate data-model residuals as a function of phase
        ###

        # phase and redshift tables
        phases_raisin = np.loadtxt('ymodelfiles/phases_RAISIN.txt',unpack=True)
        phases_csp = np.loadtxt('ymodelfiles/phases_CSP.txt',unpack=True)
        redshifts_raisin = np.loadtxt('ymodelfiles/redshifts_RAISIN.txt',unpack=True)
        redshifts_csp = np.loadtxt('ymodelfiles/redshifts_CSP.txt',unpack=True)
        # then for individual filters
        magadjust_f125w = np.loadtxt('ymodelfiles/delta_M_original-modified_F125W.txt',unpack=True)
        magadjust_f160w = np.loadtxt('ymodelfiles/delta_M_original-modified_F160W.txt',unpack=True)        
        magadjust_Y = np.loadtxt('ymodelfiles/delta_M_original-modified_Y_RC.txt',unpack=True)
        magadjust_y = np.loadtxt('ymodelfiles/delta_M_original-modified_Y_WIRC.txt',unpack=True)
        magadjust_J = np.loadtxt('ymodelfiles/delta_M_original-modified_J_RC1.txt',unpack=True)
        magadjust_j = np.loadtxt('ymodelfiles/delta_M_original-modified_J_RC2.txt',unpack=True)
        magadjust_H = np.loadtxt('ymodelfiles/delta_M_original-modified_H_RC.txt',unpack=True)
        magadjust_dict_csp = {'Y':magadjust_Y,'y':magadjust_y,'J':magadjust_J,'j':magadjust_j,'H':magadjust_H}
        magadjust_dict_raisin = {'J':magadjust_f125w,'H':magadjust_f160w}
        #import pdb; pdb.set_trace()
        
        for version,magadjust,phases,redshifts in zip(
                ['CSPDR3_RAISIN','PS1_RAISIN','DES_RAISIN'],
                [magadjust_dict_csp,magadjust_dict_raisin,magadjust_dict_raisin],
                [phases_csp,phases_raisin,phases_raisin],
                [redshifts_csp,redshifts_raisin,redshifts_raisin]):
            
            newdir = f"data/Photometry/{version}_modeladjust"
            if not os.path.exists(newdir):
                os.makedirs(newdir)
            files = glob.glob(f"data/Photometry/{version}/*")
            for f in files:
                varnames = None
                with open(f) as fin, open(f.replace(version,f"{version}_modeladjust"),"w") as fout:
                    for line in fin:
                        line = line.replace('\n','')
                        if line.startswith('VARLIST:'):
                            varnames = np.array(line.split())
                            print(line,file=fout)
                        elif line.startswith('SNID:'):
                            snid = line.split()[1]
                            print(line,file=fout)
                        elif line.startswith('PEAKMJD'):
                            pkmjd = float(line.split()[1])
                            print(line,file=fout)
                        elif line.startswith('REDSHIFT_HELIO'):
                            zhel = float(line.split()[1])
                            print(line,file=fout)
                        elif line.startswith('OBS:'):
                            lineparts = np.array(line.split())
                            flt = lineparts[(varnames == 'FLT') | (varnames == 'BAND')][0]

                            if flt not in ['Y','y','J','j','H']:
                                print(line,file=fout)
                            else:
                                ifluxcal = np.where(varnames == 'FLUXCAL')[0]
                                mjd = float(lineparts[varnames == 'MJD'][0])
                                singlephase = (mjd-pkmjd)/(1+zhel)
                                if singlephase > 40 or singlephase < -10:
                                    print(line,file=fout)
                                else:
                                    #print(singlephase)
                                    interp_phase = interp1d(phases,magadjust[flt],axis=0)
                                    magoff_phase = interp_phase(singlephase)
                                    interp_redshift = interp1d(redshifts,magoff_phase)
                                    magoff = interp_redshift(zhel)
                                    
                                    fluxcal_new = float(lineparts[ifluxcal])*10**(-0.4*magoff) # data becomes brighter to match the model
                                    #if 'PS1' in version: import pdb; pdb.set_trace()
                                    if fluxcal_new != fluxcal_new:
                                        import pdb; pdb.set_trace()
                                    lineparts[ifluxcal] = fluxcal_new
                                    print(" ".join(lineparts),file=fout)

                        else:
                            print(line,file=fout)

            os.chdir(f"data/Photometry/{version}_modeladjust")
            os.system(f'ls *dat > {version}_modeladjust.LIST')
            os.system(f'ls *DAT >> {version}_modeladjust.LIST')
            os.chdir("../../../")
            with open(f"data/Photometry/{version}_modeladjust/{version}_modeladjust.README","w") as fout:
                print("""DOCUMENTATION_START:
blah!
DOCUMENTATION_END:""",file=fout)

                
if __name__ == "__main__":
    #Yfits()
    #DESY()
    # comp DES15C3odz (late) vs. DES16S1agd (early)
    ybm = YbandModel()
    #ybm.main()
    ybm.adjust_data()
    #ybm.adjust_data_bayesn()
