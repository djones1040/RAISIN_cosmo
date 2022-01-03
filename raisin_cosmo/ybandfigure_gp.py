#!/usr/bin/env python
# D. Jones - 9/24/21
# make a template from RAISIN data, CSP data, and compare to SNooPy
# at some nominal stretch

import numpy as np
import pylab as plt
plt.ion()
import glob
import snana
from txtobj import txtobj
from snpy import kcorr
import cosmo
from astropy.io import fits
snpyfits = fits.open('snoopy.B18/SNooPy_test.fits')
import extinction
import palettable
from palettable.colorbrewer.qualitative import Dark2_8 as palettable_color
import matplotlib.colors as colors

M0s = {'u':-18.82, 'B':-19.11, 'V':-19.12, 'g':-19.16, 'r':-19.03,
       'i':-18.50, 'Y':-18.45, 'J':-18.44, 'H':-18.38, 'K':-18.42,
       'J_K':-18.44, 'H_K':-18.38,
       'Bs':-19.319, 'Vs':-19.246, 'Rs':-19.248, 'Is':-18.981,
       'UVM2':-19, 'UVW2':-19, 'UVW1':-19}

a = {0: {'B': -19.31, 'V': -19.264, 'u': -18.945, 'g': -19.345, 'r': -19.146, 'i': -18.529, 'Y': -18.532, 'J': -18.646, 'H': -18.47}, 1: {'B': -19.317, 'V': -19.278, 'u': -18.972, 'g': -19.349, 'r': -19.162, 'i': -18.55, 'Y': -18.547, 'J': -18.665, 'H': -18.49}, 2: {'B': -19.325, 'V': -19.277, 'u': -18.969, 'g': -19.359, 'r': -19.154, 'i': -18.555, 'Y': -18.56, 'J': -18.686, 'H': -18.499}, 3: {'B': -19.271, 'V': -19.246, 'g': -19.315, 'r': -19.134, 'i': -18.518, 'Y': -18.528, 'J': -18.638, 'H': -18.462}, 4: {'B': -19.276, 'V': -19.247, 'g': -19.311, 'r': -19.134, 'i': -18.524, 'Y': -18.529, 'J': -18.646, 'H': -18.477}, 5: {'B': -19.304, 'V': -19.27, 'g': -19.344, 'r': -19.154, 'i': -18.553, 'Y': -18.561, 'J': -18.687, 'H': -18.495}}
b = {0: {'B': -0.675, 'V': -0.727, 'u': -1.077, 'g': -0.719, 'r': -0.619, 'i': -0.541, 'Y': -0.387, 'J': -0.719, 'H': -0.456}, 1: {'B': -0.655, 'V': -0.718, 'u': -1.028, 'g': -0.71, 'r': -0.613, 'i': -0.53, 'Y': -0.378, 'J': -0.697, 'H': -0.431}, 2: {'B': -0.676, 'V': -0.732, 'u': -1.123, 'g': -0.719, 'r': -0.637, 'i': -0.51, 'Y': -0.35, 'J': -0.639, 'H': -0.416}, 3: {'B': -0.753, 'V': -0.791, 'g': -0.785, 'r': -0.678, 'i': -0.599, 'Y': -0.415, 'J': -0.743, 'H': -0.513}, 4: {'B': -0.73, 'V': -0.78, 'g': -0.782, 'r': -0.672, 'i': -0.589, 'Y': -0.409, 'J': -0.728, 'H': -0.489}, 5: {'B': -0.682, 'V': -0.751, 'g': -0.727, 'r': -0.655, 'i': -0.536, 'Y': -0.36, 'J': -0.633, 'H': -0.456}}
c = {0: {'B': 3.415, 'V': 2.161, 'u': 4.066, 'g': 2.76, 'r': 1.968, 'i': 0.705, 'Y': 0.232, 'J': -0.714, 'H': -0.192}, 1: {'B': 3.5, 'V': 2.249, 'u': 4.416, 'g': 2.782, 'r': 2.049, 'i': 0.848, 'Y': 0.32, 'J': -0.538, 'H': -0.005}, 2: {'B': 3.804, 'V': 2.422, 'u': 4.742, 'g': 3.098, 'r': 2.048, 'i': 1.378, 'Y': 0.975, 'J': 0.46, 'H': 0.637}, 3: {'B': 2.928, 'V': 1.867, 'g': 2.369, 'r': 1.728, 'i': 0.476, 'Y': 0.123, 'J': -0.827, 'H': -0.374}, 4: {'B': 3.053, 'V': 1.909, 'g': 2.363, 'r': 1.744, 'i': 0.557, 'Y': 0.146, 'J': -0.72, 'H': -0.188}, 5: {'B': 3.916, 'V': 2.46, 'g': 3.166, 'r': 2.155, 'i': 1.409, 'Y': 1.024, 'J': 0.639, 'H': 0.594}}

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def lighten_color(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

def MMax(band, st, calibration=0):
    '''Given self.dm15, return the absolute magnitude at maximum for the given
      filter [band].  The calibration paramter allows you to choose which
      fit (1-6) in Folatelli et al. (2009), table 9'''
    if band not in a[calibration]:
        raise ValueError("Error, filter %s cannot be fit with this calibration")

    delta = st - 1.0
        
    return a[calibration][band] + b[calibration][band]*delta +\
        c[calibration][band]*delta**2

def poly(x, x0, coefs):
   return coefs[0] + coefs[1]*(x-x0) + coefs[2]/2*(x-x0)**2

def breakpoly(x, xb, coefs, before=True):
   if before:
      return coefs[0] + coefs[1]*(x-xb) + (x<xb)*coefs[2]*(x-xb)**2
   else:
      return coefs[0] + coefs[1]*(x-xb) + (x>xb)*coefs[2]*(x-xb)**2

def deltaTmax(band,st):
    '''Given the current dm15, what is the time of maximum of [band]
      relative to B-band.'''
    if band == 'V':  return 1.33
    if band == 'u':  return -1.60
    if band == 'g':  return 0.18
    if band == 'r':  return poly(st, 1.0, [1.56, 4.24, 25.62])
    if band == 'i':  return breakpoly(st, 1.00, [-3.28, 2.02, 16.69], True)
    if band == 'Y':  return breakpoly(st, 0.95, [-4.69, -0.08, 25.43], True)
    if band == 'J':  return breakpoly(st, 0.92, [-4.03, -2.42, 14.40], True)
    if band == 'H':  return breakpoly(st, 1.02, [-4.30, 4.26, 20.40], True)
    return 0

class ybfig:
    def __init__(self):
        self.goodcids = np.concatenate((
            np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
            np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
            np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)))

    def main(self):
        # get the k-corrs for each SN at each
        # redshift that have a band that translates to y
        abs_mag_list = np.array([])
        abs_magerr_list = np.array([])
        phase_list = np.array([])
        survey_list = np.array([])
        ebv_list = np.array([])
        st_list = np.array([])
        snid_list = np.array([])
        z_list,dlmag_list,dlmagerr_list,app_mag_list,app_magerr_list = \
            np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
        
        filtdict_raisin = {'J':'f125w',
                           'H':'f160w'}
        filtdict_csp = {'Y':'Y',
                        'y':'Ydw',
                        'J':'J',
                        'j':'Jrc2',
                        'H':'H',
                        'h':'Hdw'}
        
        for lcm,fra,files,filtdict in zip(
                [txtobj(f"yband/PS1_RAISIN_optnir.LCPLOT.TEXT",fitresheader=True,rowprfx='OBS'),
                 txtobj(f"yband/DES_RAISIN_optnir.LCPLOT.TEXT",fitresheader=True,rowprfx='OBS'),
                 txtobj(f"yband/CSPDR3_RAISIN_optnir.LCPLOT.TEXT",fitresheader=True,rowprfx='OBS')],
                [txtobj(f"yband/PS1_RAISIN_optnir.FITRES.TEXT",fitresheader=True),
                 txtobj(f"yband/DES_RAISIN_optnir.FITRES.TEXT",fitresheader=True),
                 txtobj(f"yband/CSPDR3_RAISIN_optnir.FITRES.TEXT",fitresheader=True)],
                [glob.glob('data/Photometry/PS1_RAISIN/*.dat'),
                 glob.glob('data/Photometry/DES_RAISIN/*.dat'),
                 glob.glob('data/Photometry/CSPDR3_RAISIN/*.DAT')],
                [filtdict_raisin,filtdict_raisin,filtdict_csp]):
        
            for f in files: #[3:4]: #[2:3]:
                #if '2007A' not in f: continue
                sn = snana.SuperNova(f)
                if sn.SNID not in self.goodcids: continue
                z = float(sn.REDSHIFT_HELIO.split()[0])
                #fra.PKMJD[fra.CID == sn.SNID] = 56481.62
                try: phase = (sn.MJD-fra.PKMJD[fra.CID == sn.SNID][0])/(1+z)
                except: continue
                ##phase = (sn.MJD-54562.956411205137)/(1+z)
                mwebv = float(sn.MWEBV.split()[0])
                try:
                    frs = txtobj(glob.glob(f'yband/*_RAISIN_optnir_single_{sn.SNID}.FITRES.TEXT')[0],fitresheader=True)
                except:
                    continue #import pdb; pdb.set_trace()
                try:
                    fry = txtobj(glob.glob(f'yband/*_RAISIN_Y_single_{sn.SNID}.FITRES.TEXT')[0],fitresheader=True)
                except: continue
                #frs = txtobj('OUT_TEMP_40304.FITRES.TEXT',fitresheader=True)
                #fry = txtobj('OUT_TEMP_40304.FITRES.TEXT',fitresheader=True)
                #import pdb; pdb.set_trace()
                if fry.STRETCH[0] > 1.3: continue
                ebvhost = frs.AV/1.518
                #ebvhost = -0.032
                #fry.STRETCH[0] = 1.0
                #ebvhost = 0.0
                #fry.DLMAG[0] = cosmo.mu(fry.zHD)
                
                bands_rest = lcm.BAND_REST[lcm.CID == sn.SNID]
                iRest = (bands_rest == 'Y')
                obs_band = np.unique(lcm.BAND[lcm.CID == sn.SNID][iRest])
                for o in obs_band:
                    iRest = (bands_rest == 'Y') & (lcm.BAND[lcm.CID == sn.SNID] == o)
                    phase_min = np.min(lcm.Tobs[lcm.CID == sn.SNID][iRest])/(1+z)-0.01
                    phase_max = np.max(lcm.Tobs[lcm.CID == sn.SNID][iRest])/(1+z)+0.01
                    kcorr_sn = lcm.KCOR[lcm.CID == sn.SNID][iRest]
                    kcorr_sn_mask = np.ones(len(kcorr_sn))
                    #kcorr_sn,kcorr_sn_mask = kcorr.kcorr(phase[sn.FLT == o]*(1+z),'Y',filtdict[o],z,mwebv,ebvhost,R_host=1.518)
                    r_obs = kcorr.R_obs(filtdict[o],z,phase[(sn.FLT == o) & (phase <= phase_max) & (phase >= phase_min)],ebvhost,mwebv,Rv_host=1.518).flatten()
                    #r_obs = kcorr.R_obs(filtdict[o],z,phase[(sn.FLT == o)]*(1+z),ebvhost,mwebv,Rv_host=1.518,redlaw='ccm',version='H3').flatten()
                    
                    kcorr_sn_mask = np.array(kcorr_sn_mask)
                    kcorr_sn = np.array(kcorr_sn)

                    
                    mag = -2.5*np.log10(sn.FLUXCAL[(sn.FLT == o) & (phase <= phase_max) & (phase >= phase_min)][kcorr_sn_mask == 1]) + 27.5
                    magerr = 1.086*sn.FLUXCALERR[(sn.FLT == o) & (phase <= phase_max) & (phase >= phase_min)][kcorr_sn_mask == 1]\
                        /sn.FLUXCAL[(sn.FLT == o) & (phase <= phase_max) & (phase >= phase_min)][kcorr_sn_mask == 1]
                    #mag = -2.5*np.log10(sn.FLUXCAL[(sn.FLT == o)][kcorr_sn_mask == 1]) + 27.5
                    #magerr = 1.086*sn.FLUXCALERR[(sn.FLT == o)][kcorr_sn_mask == 1]\
                    #    /sn.FLUXCAL[(sn.FLT == o)][kcorr_sn_mask == 1]
                    if 'CSP' in f:
                        ebvcorr = extinction.fitzpatrick99(np.array([10350.8]),mwebv*3.1)
                        mag -= ebvcorr

                    # correct for predicted diff between template with measured stretch and template w/ mean stretch
                    try:
                        mags_st_new = self.get_snoopy(phase[(sn.FLT == o)][kcorr_sn_mask == 1],sbv=fry.STRETCH[0])
                        mags_st_1 = self.get_snoopy(phase[(sn.FLT == o)][kcorr_sn_mask == 1],sbv=1)
                        magoff = mags_st_new-mags_st_1
                    except: continue
                        
                    #magoff = np.zeros(len(mags_st_1))                    
                        
                    abs_mag_list = np.append(abs_mag_list,mag-kcorr_sn[kcorr_sn_mask == 1]-fry.DLMAG[0]-magoff-r_obs[kcorr_sn_mask == 1]*ebvhost) #cosmo.mu(z))
                    abs_magerr_list = np.append(abs_magerr_list,magerr)
                    phase_list = np.append(phase_list,phase[(sn.FLT == o) & (phase <= phase_max) & (phase >= phase_min)][kcorr_sn_mask == 1])
                    #phase_list = np.append(phase_list,phase[(sn.FLT == o)][kcorr_sn_mask == 1])
                    if 'CSP' in f:
                        survey_list = np.append(survey_list,['CSP']*len(magerr))
                    elif 'PS1' in f:
                        survey_list = np.append(survey_list,['PS1']*len(magerr))
                    elif 'DES' in f:
                        survey_list = np.append(survey_list,['DES']*len(magerr))
                    ebv_list = np.append(ebv_list,[ebvhost]*len(magerr))
                    st_list = np.append(st_list,[ebvhost]*len(magerr))
                    snid_list = np.append(snid_list,[sn.SNID]*len(magerr))
                    z_list = np.append(z_list,[fry.zCMB[0]]*len(magerr))
                    dlmag_list = np.append(dlmag_list,[fry.DLMAG[0]]*len(magerr))
                    dlmagerr_list = np.append(dlmagerr_list,[fry.DLMAGERR[0]]*len(magerr))
                    app_mag_list = np.append(app_mag_list,mag-kcorr_sn[kcorr_sn_mask == 1]-magoff-r_obs[kcorr_sn_mask == 1]*ebvhost)
                    app_magerr_list = np.append(app_magerr_list,magerr)
                    
        # plot up the data in rest-frame
        #if len(phase_list):
        with open('restmags_noebvstretch.txt','w') as fout:
            print('# SNID phase Mabs Mabserr ebvhost sbv',file=fout)
            for ebv,st,absmag,absmagerr,phase,snid in zip(ebv_list,st_list,abs_mag_list,abs_magerr_list,phase_list,snid_list):
                print(f"{snid} {phase:.3f} {absmag:.3f} {absmagerr:.3f} {ebv:.3f} {st:.3f}",file=fout)
        with open('restmags_appmag_noebvstretch.txt','w') as fout:
            print('# SN z phase m_Y m_Y_err mu_snpy mu_snpy_err mu_lcdm',file=fout)
            for ebv,st,appmag,appmagerr,phase,snid,z,dlmag,dlmagerr in zip(ebv_list,st_list,app_mag_list,app_magerr_list,phase_list,snid_list,z_list,dlmag_list,dlmagerr_list):
                print(f"{snid} {z:.5f} {phase:.3f} {appmag:.3f} {appmagerr:.3f} {dlmag:.3f} {dlmagerr:.3f} {cosmo.mu(z):.3f}",file=fout)

        def lighten_color(arg1,arg2):
            return arg1
                
        ax = plt.axes()
        ax.set_prop_cycle('color', palettable_color.mpl_colors)
        ax.errorbar(
            phase_list[survey_list == 'CSP'],
            abs_mag_list[survey_list == 'CSP'],
            yerr=abs_magerr_list[survey_list == 'CSP'],
            fmt='.',color=lighten_color('C0',0.8),alpha=0.2)
        ax.errorbar(
            phase_list[survey_list != 'CSP'],
            abs_mag_list[survey_list != 'CSP'],
            yerr=abs_magerr_list[survey_list != 'CSP'],
            fmt='.',color=lighten_color('C1',0.8),alpha=0.2)
        ax.invert_yaxis()
        ax.set_xlabel('Phase (days)',fontsize=15)
        ax.set_ylabel('$M_Y$ (mag)',fontsize=15)
        self.plot_snoopy(ax) #,sbv=1.26)#,sbv=fry.STRETCH[0])
        ax.set_ylim([-15,-19.5])

        #if np.min(abs_mag_list) < -18.7: import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()
        # do a gaussian process fit
        N = 30 #int(len(datamag[survey == version])/15) # arbitrary
        phase_sort = np.sort(phase_list[survey_list == 'CSP'])
        residmag = abs_mag_list[survey_list == 'CSP']
        residmag_sort = residmag[np.argsort(phase_list[survey_list == 'CSP'])]
        rolling_mean = np.convolve(residmag_sort, np.ones((N,))/N, mode='same')
        ax.plot(phase_sort[(phase_sort > -8) & (phase_sort < 40)],
                np.convolve(residmag_sort, np.ones((N,))/N, mode='same')[(phase_sort > -8) & (phase_sort < 40)],
                zorder=100,lw=2,label='CSP',color=lighten_color('C0',1.2))

        N = 10 #int(len(datamag[survey == version])/15) # arbitrary
        phase_sort = np.sort(phase_list[survey_list != 'CSP'])
        residmag = abs_mag_list[survey_list != 'CSP']
        residmag_sort = residmag[np.argsort(phase_list[survey_list != 'CSP'])]
        rolling_mean = np.convolve(residmag_sort, np.ones((N,))/N, mode='same')
        ax.plot(phase_sort[(phase_sort > 6) & (phase_sort < 25)],
                np.convolve(residmag_sort, np.ones((N,))/N, mode='same')[(phase_sort > 6) & (phase_sort < 25)],
                zorder=100,lw=2,label='RAISIN',color=lighten_color('C1',1.2))

        ax.set_xlim([-10,40])
        ax.set_ylim([-17.25,-18.75])
        ax.legend(loc='upper center')
        
        # plot the fit

        import pdb; pdb.set_trace()

    def main_kaiseyfig(self):
        # get the k-corrs for each SN at each
        # redshift that have a band that translates to y
        abs_mag_list = np.array([])
        abs_magerr_list = np.array([])
        phase_list = np.array([])
        survey_list = np.array([])
        ebv_list = np.array([])
        st_list = np.array([])
        snid_list = np.array([])
        z_list,dlmag_list,dlmagerr_list,app_mag_list,app_magerr_list = \
            np.array([]),np.array([]),np.array([]),np.array([]),np.array([])

        from matplotlib.gridspec import GridSpec
        fig = plt.figure()
        gs = GridSpec(4,6, figure=fig)
        gs.update(wspace=1.0,hspace=0.0) #6)
        ax1 = fig.add_subplot(gs[0:3,0:5])
        ax2 = fig.add_subplot(gs[3,0:5])
        cax = fig.add_subplot(gs[:,5])
        for ax in [ax1,ax2,cax]:
            ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        cax.yaxis.tick_right()
        #import pdb; pdb.set_trace()
        filtdict_raisin = {'J':'f125w',
                           'H':'f160w'}
        filtdict_csp = {'Y':'Y',
                        'y':'Ydw',
                        'J':'J',
                        'j':'Jrc2',
                        'H':'H',
                        'h':'Hdw'}
        
        for lcm,fra,files,filtdict in zip(
                [txtobj(f"yband/PS1_RAISIN_optnir.LCPLOT.TEXT",fitresheader=True,rowprfx='OBS'),
                 txtobj(f"yband/DES_RAISIN_optnir.LCPLOT.TEXT",fitresheader=True,rowprfx='OBS'),
                 txtobj(f"yband/CSPDR3_RAISIN_optnir.LCPLOT.TEXT",fitresheader=True,rowprfx='OBS')],
                [txtobj(f"yband/PS1_RAISIN_optnir.FITRES.TEXT",fitresheader=True),
                 txtobj(f"yband/DES_RAISIN_optnir.FITRES.TEXT",fitresheader=True),
                 txtobj(f"yband/CSPDR3_RAISIN_optnir.FITRES.TEXT",fitresheader=True)],
                [glob.glob('data/Photometry/PS1_RAISIN/*.dat'),
                 glob.glob('data/Photometry/DES_RAISIN/*.dat'),
                 glob.glob('data/Photometry/CSPDR3_RAISIN/*.DAT')],
                [filtdict_raisin,filtdict_raisin,filtdict_csp]):
        
            for f in files: #[3:4]: #[2:3]:
                #if '2007A' not in f: continue
                sn = snana.SuperNova(f)
                if sn.SNID not in self.goodcids: continue
                z = float(sn.REDSHIFT_HELIO.split()[0])
                #fra.PKMJD[fra.CID == sn.SNID] = 56481.62
                try: phase = (sn.MJD-fra.PKMJD[fra.CID == sn.SNID][0])/(1+z)
                except: continue
                ##phase = (sn.MJD-54562.956411205137)/(1+z)
                mwebv = float(sn.MWEBV.split()[0])
                try:
                    frs = txtobj(glob.glob(f'yband/*_RAISIN_optnir_single_{sn.SNID}.FITRES.TEXT')[0],fitresheader=True)
                except:
                    continue #import pdb; pdb.set_trace()
                try:
                    fry = txtobj(glob.glob(f'yband/*_RAISIN_Y_single_{sn.SNID}.FITRES.TEXT')[0],fitresheader=True)
                except: continue
                #frs = txtobj('OUT_TEMP_40304.FITRES.TEXT',fitresheader=True)
                #fry = txtobj('OUT_TEMP_40304.FITRES.TEXT',fitresheader=True)
                #import pdb; pdb.set_trace()
                if fry.STRETCH[0] > 1.3: continue
                ebvhost = frs.AV/1.518
                #ebvhost = -0.032
                fry.STRETCH[0] = 1.0
                ebvhost = 0.0
                #fry.DLMAG[0] = cosmo.mu(fry.zHD)
                
                bands_rest = lcm.BAND_REST[lcm.CID == sn.SNID]
                iRest = (bands_rest == 'Y')
                obs_band = np.unique(lcm.BAND[lcm.CID == sn.SNID][iRest])
                for o in obs_band:
                    iRest = (bands_rest == 'Y') & (lcm.BAND[lcm.CID == sn.SNID] == o)
                    phase_min = np.min(lcm.Tobs[lcm.CID == sn.SNID][iRest])/(1+z)-0.01
                    phase_max = np.max(lcm.Tobs[lcm.CID == sn.SNID][iRest])/(1+z)+0.01
                    kcorr_sn = lcm.KCOR[lcm.CID == sn.SNID][iRest]
                    kcorr_sn_mask = np.ones(len(kcorr_sn))
                    #kcorr_sn,kcorr_sn_mask = kcorr.kcorr(phase[sn.FLT == o]*(1+z),'Y',filtdict[o],z,mwebv,ebvhost,R_host=1.518)
                    r_obs = kcorr.R_obs(filtdict[o],z,phase[(sn.FLT == o) & (phase <= phase_max) & (phase >= phase_min)],ebvhost,mwebv,Rv_host=1.518).flatten()
                    #r_obs = kcorr.R_obs(filtdict[o],z,phase[(sn.FLT == o)]*(1+z),ebvhost,mwebv,Rv_host=1.518,redlaw='ccm',version='H3').flatten()
                    
                    kcorr_sn_mask = np.array(kcorr_sn_mask)
                    kcorr_sn = np.array(kcorr_sn)

                    
                    mag = -2.5*np.log10(sn.FLUXCAL[(sn.FLT == o) & (phase <= phase_max) & (phase >= phase_min)][kcorr_sn_mask == 1]) + 27.5
                    magerr = 1.086*sn.FLUXCALERR[(sn.FLT == o) & (phase <= phase_max) & (phase >= phase_min)][kcorr_sn_mask == 1]\
                        /sn.FLUXCAL[(sn.FLT == o) & (phase <= phase_max) & (phase >= phase_min)][kcorr_sn_mask == 1]
                    #mag = -2.5*np.log10(sn.FLUXCAL[(sn.FLT == o)][kcorr_sn_mask == 1]) + 27.5
                    #magerr = 1.086*sn.FLUXCALERR[(sn.FLT == o)][kcorr_sn_mask == 1]\
                    #    /sn.FLUXCAL[(sn.FLT == o)][kcorr_sn_mask == 1]
                    if 'CSP' in f:
                        ebvcorr = extinction.fitzpatrick99(np.array([10350.8]),mwebv*3.1)
                        mag -= ebvcorr

                    # correct for predicted diff between template with measured stretch and template w/ mean stretch
                    try:
                        mags_st_new = self.get_snoopy(phase[(sn.FLT == o)][kcorr_sn_mask == 1],sbv=fry.STRETCH[0])
                        mags_st_1 = self.get_snoopy(phase[(sn.FLT == o)][kcorr_sn_mask == 1],sbv=1)
                        magoff = mags_st_new-mags_st_1
                    except: continue
                        
                    magoff = np.zeros(len(mags_st_1))                    
                        
                    abs_mag_list = np.append(abs_mag_list,mag-kcorr_sn[kcorr_sn_mask == 1]-cosmo.mu(fry.zHD[0])-magoff-r_obs[kcorr_sn_mask == 1]*ebvhost) #fry.DLMAG[0]#cosmo.mu(z))
                    abs_magerr_list = np.append(abs_magerr_list,magerr)
                    phase_list = np.append(phase_list,phase[(sn.FLT == o) & (phase <= phase_max) & (phase >= phase_min)][kcorr_sn_mask == 1])
                    #phase_list = np.append(phase_list,phase[(sn.FLT == o)][kcorr_sn_mask == 1])
                    if 'CSP' in f:
                        survey_list = np.append(survey_list,['CSP']*len(magerr))
                    elif 'PS1' in f:
                        survey_list = np.append(survey_list,['PS1']*len(magerr))
                    elif 'DES' in f:
                        survey_list = np.append(survey_list,['DES']*len(magerr))
                    ebv_list = np.append(ebv_list,[ebvhost]*len(magerr))
                    st_list = np.append(st_list,[ebvhost]*len(magerr))
                    snid_list = np.append(snid_list,[sn.SNID]*len(magerr))
                    z_list = np.append(z_list,[fry.zCMB[0]]*len(magerr))
                    dlmag_list = np.append(dlmag_list,[fry.DLMAG[0]]*len(magerr))
                    dlmagerr_list = np.append(dlmagerr_list,[fry.DLMAGERR[0]]*len(magerr))
                    app_mag_list = np.append(app_mag_list,mag-kcorr_sn[kcorr_sn_mask == 1]-magoff-r_obs[kcorr_sn_mask == 1]*ebvhost)
                    app_magerr_list = np.append(app_magerr_list,magerr)
                    
        # plot up the data in rest-frame
        #if len(phase_list):
        #with open('restmags_noebvstretch.txt','w') as fout:
        #    print('# SNID phase Mabs Mabserr ebvhost sbv',file=fout)
        #    for ebv,st,absmag,absmagerr,phase,snid in zip(ebv_list,st_list,abs_mag_list,abs_magerr_list,phase_list,snid_list):
        #        print(f"{snid} {phase:.3f} {absmag:.3f} {absmagerr:.3f} {ebv:.3f} {st:.3f}",file=fout)
        #with open('restmags_appmag_noebvstretch.txt','w') as fout:
        #    print('# SN z phase m_Y m_Y_err mu_snpy mu_snpy_err mu_lcdm',file=fout)
        #    for ebv,st,appmag,appmagerr,phase,snid,z,dlmag,dlmagerr in zip(ebv_list,st_list,app_mag_list,app_magerr_list,phase_list,snid_list,z_list,dlmag_list,dlmagerr_list):
        #        print(f"{snid} {z:.5f} {phase:.3f} {appmag:.3f} {appmagerr:.3f} {dlmag:.3f} {dlmagerr:.3f} {cosmo.mu(z):.3f}",file=fout)
        
        #ax = plt.axes()
        ax1.set_prop_cycle('color', palettable_color.mpl_colors)
        #ax1.errorbar(
        #    phase_list[survey_list == 'CSP'],
        #    abs_mag_list[survey_list == 'CSP'],
        #    yerr=abs_magerr_list[survey_list == 'CSP'],
        #    fmt='.',color=lighten_color('C0',0.8), alpha=0.5)
        import matplotlib.cm as cm
        from matplotlib.colors import Normalize
        #cmap = cm.winter
        import palettable
        cmap = truncate_colormap(cm.viridis,minval=0.2,maxval=1.0) #palettable.scientific.diverging.Berlin_15.mpl_colormap #palettable.colorbrewer.qualitative.Dark2_4.mpl_colormap #palettable.cartocolors.sequential.DarkMint_7.mpl_colormap #palettable.cartocolors.sequential.agSunset_7.mpl_colormap #palettable_color.mpl_colormap #cm.brg #terrain #palettable.scientific.sequential.Imola_16.mpl_colormap
        def lighten_color(arg1,arg2):
            return arg1

        norm = Normalize(vmin=z_list[survey_list != 'CSP'].min(), vmax=z_list[survey_list != 'CSP'].max())
        ax1.errorbar(
            phase_list, #[survey_list != 'CSP'],
            abs_mag_list, #[survey_list != 'CSP'],
            yerr=abs_magerr_list, #[survey_list != 'CSP'],
            fmt='none', ecolor=cmap(norm(z_list)),alpha=0.8) #[survey_list != 'CSP'])))
        sc = ax1.scatter(
            phase_list, #[survey_list != 'CSP'],
            abs_mag_list, #[survey_list != 'CSP'],
            c=z_list, #[survey_list != 'CSP'],
            s=10,marker='o',
            cmap=cmap,alpha=0.8)
        ax1.invert_yaxis()
        #ax1.set_xlabel('Phase (days)',fontsize=15)
        ax1.set_ylabel('$m_Y - \mu_{\Lambda CDM}$ (mag)',fontsize=15)
        self.plot_snoopy(ax1) #,sbv=1.26)#,sbv=fry.STRETCH[0])
        ax1.set_ylim([-15,-19.5])

        #if np.min(abs_mag_list) < -18.7: import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()
        # do a gaussian process fit
        N = 30 #int(len(datamag[survey == version])/15) # arbitrary
        phase_sort = np.sort(phase_list[survey_list == 'CSP'])
        zlist = z_list[survey_list == 'CSP']
        z_sort = z_list[survey_list == 'CSP'][np.argsort(phase_list[survey_list == 'CSP'])]
        residmagerr_sort = abs_magerr_list[survey_list == 'CSP'][np.argsort(phase_list[survey_list == 'CSP'])]
        residmag = abs_mag_list[survey_list == 'CSP']
        residmagerr = abs_magerr_list[survey_list == 'CSP']
        residmag_sort = residmag[np.argsort(phase_list[survey_list == 'CSP'])]
        rolling_mean = np.convolve(residmag_sort, np.ones((N,))/N, mode='same')
        ax1.plot(phase_sort[(phase_sort > -8) & (phase_sort < 40)],
                np.convolve(residmag_sort, np.ones((N,))/N, mode='same')[(phase_sort > -8) & (phase_sort < 40)],
                zorder=100,lw=2,label='CSP',color=lighten_color('C0',0.8)) #palettable_color.mpl_colors[0]) #lighten_color('C0',1.2))

        N = 10 #int(len(datamag[survey == version])/15) # arbitrary
        phase_sort = np.sort(phase_list[survey_list != 'CSP'])
        residmag = abs_mag_list[survey_list != 'CSP']
        residmagerr = abs_magerr_list[survey_list != 'CSP']
        zlist = z_list[survey_list != 'CSP']
        residmag_sort = residmag[np.argsort(phase_list[survey_list != 'CSP'])]
        rolling_mean = np.convolve(residmag_sort, np.ones((N,))/N, mode='same')
        ax1.plot(phase_sort[(phase_sort > 6) & (phase_sort < 25)],
                np.convolve(residmag_sort, np.ones((N,))/N, mode='same')[(phase_sort > 6) & (phase_sort < 25)],
                zorder=100,lw=2,label='RAISIN',color=lighten_color(cmap(0.7),1.2)) #palettable_color.mpl_colors[1]) #lighten_color('C1',1.2))

        ax1.set_xlim([-10,40])
        ax1.set_ylim([-16.9,-19])
        ax1.legend(loc='upper center')
        
        # plot the variance
        import pandas as pd
        phase_sort = np.sort(phase_list)
        residmag = abs_mag_list[:]
        residmagerr = abs_magerr_list[:]
        zlist = z_list[:]
        residmag_sort = residmag[np.argsort(phase_list)]
        residmagerr_sort = residmagerr[np.argsort(phase_list)]
        z_sort = zlist[np.argsort(phase_list)]
        
        ts = pd.Series(residmag_sort)
        rolling_std = ts.rolling(window=50).std()
        def rolling_std(z_sort,phase_sort,residmag_sort,residmagerr_sort,window=50):
            std_adj = np.array([])
            for i,r in enumerate(residmag_sort):
                ilow = i-window if i-window > 0 else 0
                ihi = i+window if i+window < len(residmag_sort) else len(residmag_sort)
                std = np.std(residmag_sort[ilow:ihi])
                z = np.median(z_sort[ilow:ihi])
                peczerr = 250/3e5
                zerr = peczerr*5.0/np.log(10)*(1.0+z)/(z*(1.0+z/2.0))
                err = np.median(residmagerr_sort[ilow:ihi])
                std_adj = np.append(std_adj,np.sqrt(std**2.-zerr**2.-err**2.))
                #import pdb; pdb.set_trace()
            return std_adj

        iphase = [(phase_sort > -4) & (phase_sort < 40)]
        #import pdb; pdb.set_trace()
        rstd = rolling_std(z_sort[iphase],phase_sort[iphase],residmag_sort[iphase],residmagerr_sort[iphase],window=25)
        ax2.plot(phase_sort[(phase_sort > -4) & (phase_sort < 40)],rstd,
                 color=lighten_color('C0',0.8)) #palettable_color.mpl_colors[0])

        phase_sort = np.sort(phase_list[survey_list != 'CSP'])
        residmag = abs_mag_list[survey_list != 'CSP']
        residmag_sort = residmag[np.argsort(phase_list[survey_list != 'CSP'])]
        residmagerr = abs_magerr_list[survey_list != 'CSP']
        zlist = z_list[survey_list != 'CSP']
        residmag_sort = residmag[np.argsort(phase_list[survey_list != 'CSP'])]
        residmagerr_sort = residmagerr[np.argsort(phase_list[survey_list != 'CSP'])]
        z_sort = zlist[np.argsort(phase_list[survey_list != 'CSP'])]

        N = 10
        ts = pd.Series(residmag_sort)
        #rolling_std = ts.rolling(window=30).std()

        iphase = [(phase_sort > 6) & (phase_sort < 25)]
        rstd = rolling_std(z_sort[iphase],phase_sort[iphase],residmag_sort[iphase],residmagerr_sort[iphase],window=15)
        ax2.plot(phase_sort[(phase_sort > 6) & (phase_sort < 25)],rstd, #rolling_std[(phase_sort > 6) & (phase_sort < 25)],
                 color=lighten_color(cmap(0.7),1.2)) #palettable_color.mpl_colors[1])

        
        ax2.set_xlim([-10,40])
        ax2.set_xlabel('Phase (days)',fontsize=15)
        ax2.set_ylabel('$\sigma_{int}$ (mag)',fontsize=15)
        ax2.set_ylim([0.1,0.29])

        
        cb = plt.colorbar(sc,cax=cax,orientation='vertical')
        cb.set_alpha(1.0)
        cb.draw_all()
        cax.set_ylabel('$z_{CMB}$',fontsize=15)
        ax1.yaxis.set_ticks([-17.0,-17.5,-18.0,-18.5,-19.0])
        ax1.xaxis.set_ticklabels([])
        
        import pdb; pdb.set_trace()

        
    def plot_snoopy(self,ax,sbv=1):

        filters = snpyfits[0].header['FILTERS']

        NGRID_FILT = snpyfits[1].data['NBIN'][4]
        NGRID_TREST = snpyfits[1].data['NBIN'][5]
        NWDPAD = 4
        ST_GRID = np.linspace(0.7,1.3,snpyfits[1].data['NBIN'][snpyfits[1].data['PAR_NAME'] == 'STRETCH'][0])

        I2MAG = snpyfits[9].data['I2MAG']
        I2MAGERR = snpyfits[9].data['I2MAGERR']
        TREST = snpyfits[7].data['TREST']

        ILCOFF = [1, 1, 1, 1]
        #norm=colors.Normalize(vmin=0.7,vmax=1.3)
        cmap=plt.get_cmap('RdBu')
        for st in ST_GRID:

            # redshift, AV, RV, luminosity
            INDX = [1,1,1,np.where(ST_GRID == st)[0]+1]
            ILC = 1
            for i in range(4):
                ILC += ILCOFF[i] * (INDX[i] - 1)

            PTR = 1 + ((NGRID_FILT * NGRID_TREST) + NWDPAD) * (ILC-1)

            for iflt,flt in enumerate(filters):

                lctmpl = snpyfits[9].data['I2MAG'][int(PTR)+1+NGRID_TREST*iflt:int(PTR)+1+NGRID_TREST*(iflt+1)]/1000.
                #color=cmap(norm(st))
                if flt == 'Y' and st > sbv-0.002 and st < sbv+0.002: ax.plot(TREST,lctmpl,color='k',lw=2,label='SNooPy Template')

    def get_snoopy(self,phase,sbv=1):

        snpyfits = fits.open('snoopy.B18/SNooPy_test.fits')
        filters = snpyfits[0].header['FILTERS']

        NGRID_FILT = snpyfits[1].data['NBIN'][4]
        NGRID_TREST = snpyfits[1].data['NBIN'][5]
        NWDPAD = 4
        ST_GRID = np.linspace(0.7,1.3,snpyfits[1].data['NBIN'][snpyfits[1].data['PAR_NAME'] == 'STRETCH'][0])

        I2MAG = snpyfits[9].data['I2MAG']
        I2MAGERR = snpyfits[9].data['I2MAGERR']
        TREST = snpyfits[7].data['TREST']

        ILCOFF = [1, 1, 1, 1]
        #norm=colors.Normalize(vmin=0.7,vmax=1.3)
        cmap=plt.get_cmap('RdBu')
        for st in ST_GRID:

            # redshift, AV, RV, luminosity
            INDX = [1,1,1,np.where(ST_GRID == st)[0]+1]
            ILC = 1
            for i in range(4):
                ILC += ILCOFF[i] * (INDX[i] - 1)

            PTR = 1 + ((NGRID_FILT * NGRID_TREST) + NWDPAD) * (ILC-1)

            for iflt,flt in enumerate(filters):

                lctmpl = snpyfits[9].data['I2MAG'][int(PTR)+1+NGRID_TREST*iflt:int(PTR)+1+NGRID_TREST*(iflt+1)]/1000.
                #color=cmap(norm(st))
                if flt == 'Y' and st > sbv-0.002 and st < sbv+0.002:
                    return np.interp(phase,TREST,lctmpl)
                    
        
if __name__ == "__main__":
    yb = ybfig()
    #yb.main()
    yb.main_kaiseyfig()
