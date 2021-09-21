#!/usr/bin/env python
# D. Jones - 9/20/21

import numpy as np
from txtobj import txtobj
import pylab as plt
plt.ion()
from scipy.optimize import minimize
import cosmo
import glob
import snana

import palettable
from palettable.colorbrewer.qualitative import Dark2_8 as palettable_color

def main():

    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)

    # plot the residuals
    plot_resids(ax1)
    
    # plot the y-band model diff
    yp = YbandModelPlot()
    yp.ymodeldiff(ax2)

    ax1.set_xlabel('Median NIR Phase')
    ax2.set_xlabel('Phase')

    ax1.set_ylabel('High-$z$ Hubble Resid. (mag)')
    ax2.set_ylabel('Rest-Frame $Y$ Residual (mag)')
    
    import pdb; pdb.set_trace()
    
def plot_resids(ax):

    snidlist = []
    min_phase_list = []
    med_phase_list = []
    max_phase_list = []
    host_logmass_list = []
    nfilt_list = []
    files = glob.glob('data/Photometry/PS1_RAISIN/*.snana.dat')
    for f in files:
        sn = snana.SuperNova(f)
        min_phase_list += [(np.min(sn.MJD[(sn.FLT == 'J') | (sn.FLT == 'H')])-sn.PEAKMJD)/(1+float(sn.REDSHIFT_HELIO.split()[0]))]
        med_phase_list += [(np.median(sn.MJD[(sn.FLT == 'J') | (sn.FLT == 'H')])-sn.PEAKMJD)/(1+float(sn.REDSHIFT_HELIO.split()[0]))]
        max_phase_list += [(np.max(sn.MJD[(sn.FLT == 'J') | (sn.FLT == 'H')])-sn.PEAKMJD)/(1+float(sn.REDSHIFT_HELIO.split()[0]))]
        nfilt_list += [len(np.unique(sn.FLT[(sn.FLT == 'J') | (sn.FLT == 'H')]))]
        host_logmass_list += [float(sn.HOSTGAL_LOGMASS.split()[0])]
        snidlist += [sn.SNID]
        
    files2 = glob.glob('data/Photometry/DES_RAISIN/*.snana.dat')
    for f in files2:
        sn = snana.SuperNova(f)
        min_phase_list += [(np.min(sn.MJD[(sn.FLT == 'J') | (sn.FLT == 'H')])-sn.PEAKMJD)/(1+float(sn.REDSHIFT_HELIO.split()[0]))]
        med_phase_list += [(np.median(sn.MJD[(sn.FLT == 'J') | (sn.FLT == 'H')])-sn.PEAKMJD)/(1+float(sn.REDSHIFT_HELIO.split()[0]))]
        max_phase_list += [(np.max(sn.MJD[(sn.FLT == 'J') | (sn.FLT == 'H')])-sn.PEAKMJD)/(1+float(sn.REDSHIFT_HELIO.split()[0]))]
        nfilt_list += [len(np.unique(sn.FLT[(sn.FLT == 'J') | (sn.FLT == 'H')]))]
        host_logmass_list += [float(sn.HOSTGAL_LOGMASS.split()[0])]
        snidlist += [sn.SNID]
        
    min_phase_list = np.array(min_phase_list)
    med_phase_list = np.array(med_phase_list)
    max_phase_list = np.array(max_phase_list)
    host_logmass_list = np.array(host_logmass_list)
    nfilt_list = np.array(nfilt_list)
    snidlist = np.array(snidlist)
    
    frsnana = txtobj('output/fit_nir/PS1_RAISIN.FITRES.TEXT',fitresheader=True)
    frsnana2 = txtobj('output/fit_nir/DES_RAISIN.FITRES.TEXT',fitresheader=True)

    for k in frsnana.__dict__.keys():
        frsnana.__dict__[k] = np.append(frsnana.__dict__[k],frsnana2.__dict__[k])

    frsnanaopt = txtobj('output/fit_optical/PS1_RAISIN_optical.FITRES.TEXT',fitresheader=True)
    frsnanaopt2 = txtobj('output/fit_optical/DES_RAISIN_optical.FITRES.TEXT',fitresheader=True)

    for k in frsnanaopt.__dict__.keys():
        frsnanaopt.__dict__[k] = np.append(frsnanaopt.__dict__[k],frsnanaopt2.__dict__[k])

        
    goodcids = np.append(np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',dtype=str),
                         np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',dtype=str))
    iGood = np.array([],dtype=int)
    frsnana.min_phase = np.zeros(len(frsnana.CID))
    frsnana.med_phase = np.zeros(len(frsnana.CID))
    frsnana.max_phase = np.zeros(len(frsnana.CID))
    frsnana.nfilt = np.zeros(len(frsnana.CID))
    frsnana.host_logmass = np.zeros(len(frsnana.CID))
    frsnana.stretch = np.zeros(len(frsnana.CID))
    for j,i in enumerate(frsnana.CID):
        if i in goodcids:
            iGood = np.append(iGood,j)
            frsnana.min_phase[j] = min_phase_list[snidlist == i][0]
            frsnana.med_phase[j] = med_phase_list[snidlist == i][0]
            frsnana.max_phase[j] = max_phase_list[snidlist == i][0]
            frsnana.nfilt[j] = nfilt_list[snidlist == i][0]
            try: frsnana.stretch[j] = frsnanaopt.STRETCH[frsnanaopt.CID == i][0]
            except:
                import pdb; pdb.set_trace()
                #frsnana.stretch[j] = -99
            try: frsnana.host_logmass[j] =  host_logmass_list[snidlist == i][0]
            except: import pdb; pdb.set_trace()
            
    for k in frsnana.__dict__.keys():
        frsnana.__dict__[k] = frsnana.__dict__[k][iGood]

    ax.errorbar(frsnana.med_phase[frsnana.med_phase > 15],
                frsnana.DLMAG[frsnana.med_phase > 15]-cosmo.mu(frsnana.zHD[frsnana.med_phase > 15]),
                yerr=frsnana.DLMAGERR[frsnana.med_phase > 15],fmt='o',color='C0')
    ax.errorbar(frsnana.med_phase[frsnana.med_phase < 15],
                frsnana.DLMAG[frsnana.med_phase < 15]-cosmo.mu(frsnana.zHD[frsnana.med_phase < 15]),
                yerr=frsnana.DLMAGERR[frsnana.med_phase < 15],fmt='o',color='C1')
    frsnana.resid = frsnana.DLMAG-cosmo.mu(frsnana.zHD)
    frsnana.resid -= np.median(frsnana.resid)
    
    # now get the maximum likelihood stuff
    md1 = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                  args=(frsnana.resid[frsnana.med_phase < 15],frsnana.DLMAGERR[frsnana.med_phase < 15],None))
    md2 = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                  args=(frsnana.resid[frsnana.med_phase > 15],frsnana.DLMAGERR[frsnana.med_phase > 15],None))

    ax.plot([5,15],[md1.x[0],md1.x[0]],color='C1')
    ax.plot([15,25],[md2.x[0],md2.x[0]],color='C0')
    phasevals1 = np.linspace(0,15,10)
    ax.fill_between(phasevals1,[md1.x[0]-np.sqrt(md1.hess_inv[0,0])]*len(phasevals1),
                    [md1.x[0]+np.sqrt(md1.hess_inv[0,0])]*len(phasevals1),color='C1',alpha=0.2)
    phasevals2 = np.linspace(15,30,10)
    ax.fill_between(phasevals2,[md2.x[0]-np.sqrt(md2.hess_inv[0,0])]*len(phasevals2),
                    [md2.x[0]+np.sqrt(md2.hess_inv[0,0])]*len(phasevals2),color='C0',alpha=0.2)
    ax.axvline(15,color='0.8',ls='--')
    ax.axhline(0,color='k')

    ax.set_xlim(6,23)
    diff = md2.x[0]-md1.x[0]
    differr = np.sqrt(md1.hess_inv[0,0]+md2.hess_inv[0,0])
    ax.text(0.02,0.95,fr'$\Delta_{{\rm HR}} = {diff:.3f}\pm{differr:.3f}$',ha='left',va='top',transform=ax.transAxes)
    ax.set_ylim([-0.42,0.42])
    # 0.205 +/- 0.039
    # 1.100935 0.97895
    
    import pdb; pdb.set_trace()
    
    return

def lnlikefunc(x,mu_i=None,sigma_i=None,sigma=None,z=None,survey=None):

    if sigma or sigma == 0.0:
        # fix the dispersion
        x[2] = sigma; x[3] = sigma

    return -np.sum(-(mu_i-x[0])**2./(2.0*(sigma_i**2.+x[2]**2.)) +\
                np.log(1/(np.sqrt(2*np.pi)*np.sqrt(x[2]**2.+sigma_i**2.))))

class YbandModelPlot:
    def __init__(self):
        self.band = 'Y'
        self.versions = ['CSPDR3_RAISIN','PS1_RAISIN','DES_RAISIN']
        self.CIDS = [np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                     np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                     np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)]
        self.kcors = ['kcor_CSPDR3_BD17.fits','kcor_PS1MD_NIR.fits','kcor_DES_NIR.fits']
        self.filtlist = ['BomngriYyJjH','grizJH','grizJH']


    def ymodeldiff(self,ax):

        # 1. get CID list
        # 2. run opt+NIR SNANA fits
        phase,dataflux,datafluxerr,modelflux,survey,cidout = np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
        stretchout,colorout,zout = np.array([]),np.array([]),np.array([])
        for cids,kcor,filtlist,version,title in zip(self.CIDS,self.kcors,self.filtlist,self.versions,['CSP','PS1 (RAISIN1)','DES (RAISIN2)']):
            #if 'CSP' in version: continue
            outfile = f"{self.band.lower()}band/{version}_optnir"

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
                iRest = (bands_rest == self.band)
                bands_to_avoid = np.unique(lcm.BAND[lcm.CID == c][iRest])
                cid_filtdict[c] = bands_to_avoid



                # 4. run opt+NIR SNANA fits w/o that band (probably do one SN at a time)
                pardict = {}
                outfile = f"{self.band.lower()}band/{version}_optnir_single_{c}"
                fr = txtobj(f"{outfile}.FITRES.TEXT",fitresheader=True)
                try: pardict[c] = (fr.DLMAG[0],fr.STRETCH[0],fr.AV[0],fr.zHEL[0])
                except: continue
                
                outfile = f"{self.band.lower()}band/{version}_{self.band}_single_{c}_nodlmag"
                frmt = txtobj(f"{outfile}.FITRES.TEXT",fitresheader=True)
                if fr.STRETCH[0] <1.05: continue
                if fr.AV[0] > 0.4: continue

                # 5. run SNANA fits w/ Y-band only
                outfile = f"{self.band.lower()}band/{version}_{self.band}_single_{c}"

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
                    modelflux = np.append(modelflux,mf)
                    survey = np.append(survey,version)
                    cidout = np.append(cidout,c)
                    stretchout = np.append(stretchout,fr.STRETCH[0])
                    colorout = np.append(colorout,fr.AV[0])
                    zout = np.append(zout,fr.zHD[0])


            datamag = -2.5*np.log10(dataflux)+27.5
            datamagerr = 1.086*datafluxerr/dataflux
            modelmag = -2.5*np.log10(modelflux)+27.5

        # 8. run non-parametric fit to get suggested model adjustments, weighting each SN equally somehow
        #    or a rolling median?
        N = int(len(datamag)/15) # arbitrary
        phase_sort = np.sort(phase)
        residmag = datamag-modelmag
        residmag_sort = residmag[np.argsort(phase)]

        # 7. compare data to model
        datamag = -2.5*np.log10(dataflux)+27.5
        datamagerr = 1.086*datafluxerr/dataflux
        modelmag = -2.5*np.log10(modelflux)+27.5
        ax.axhline(0,color='k')
        ax.set_xlabel('Phase')
        ax.set_ylabel(f"${self.band}$ data-model")
        ax.legend()

        # 8. run non-parametric fit to get suggested model adjustments, weighting each SN equally somehow
        #    or a rolling median?
        N = int(len(datamag)/15) # arbitrary
        phase_sort = np.sort(phase)
        residmag = datamag-modelmag
        residmag_sort = residmag[np.argsort(phase)]
        survey_sort = survey[np.argsort(phase)]
        ax.set_prop_cycle('color', palettable_color.mpl_colors)
        N = 50
        ax.plot(phase_sort[survey_sort == 'CSPDR3_RAISIN'],
                    np.convolve(residmag_sort[survey_sort == 'CSPDR3_RAISIN'], np.ones((N,))/N, mode='same'),
                    zorder=100,lw=4,label='CSP',color=palettable_color.mpl_colors[0])
        N = 10
        ax.plot(phase_sort[survey_sort != 'CSPDR3_RAISIN'],
                    np.convolve(residmag_sort[survey_sort != 'CSPDR3_RAISIN'],np.ones((N,))/N, mode='same'),
                    zorder=100,lw=4,label='RAISIN',color=palettable_color.mpl_colors[1])

        # 9. now get the "final" model adjustment by taking a straight average of the low- and high-z results
        #    though we might want to test using only the high-z results to see if Kaisey's trend disappears
        phase_bins = np.arange(-20,50,2.0)


        ax.legend()
        ax.set_xlim(6,23)
        ax.set_ylim([-0.42,0.42])        
        import pdb; pdb.set_trace()
        
        return
        
if __name__ == "__main__":
    main()
