#!/usr/bin/env python
# D. Jones - 2/27/20

import numpy as np
from scipy.stats import binned_statistic
from txtobj import txtobj
import pylab as plt
plt.ion()
import glob
import snana
import os

from raisin_cosmo import ovdatamc

# rsync -avz output/sim/PS1_RAISIN_SIM_NIR david@turanga.ucsc.edu:/Users/David/Dropbox/research/RAISIN/cosmo/output/sim/

class biascor:
	def __init__(self):

		self.nirsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN_NIR_SIM',
								 '$RAISIN_ROOT/cosmo/output/fit_nir/CfA_RAISIN_NIR_SIM',
								 '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN_NIR_SIM',
								 '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN_NIR_SIM']
		self.nirsimfitreslist = [os.path.expandvars(filepath) for filepath in self.nirsimfitreslist]
		
		self.opticalsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_all/CSP_RAISIN_OPTNIR_SIM',
									 '$RAISIN_ROOT/cosmo/output/fit_all/CfA_RAISIN_OPTNIR_SIM',
									 '$RAISIN_ROOT/cosmo/output/fit_all/PS1_RAISIN_OPTNIR_SIM',
									 '$RAISIN_ROOT/cosmo/output/fit_all/DES_RAISIN_OPTNIR_SIM']
		self.opticalsimfitreslist = [os.path.expandvars(filepath) for filepath in self.opticalsimfitreslist]
		
		self.nirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN.FITRES.TEXT',
								  '$RAISIN_ROOT/cosmo/output/fit_nir/CfA_RAISIN.FITRES.TEXT',
								  '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN.FITRES.TEXT',
								  '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN.FITRES.TEXT']
		self.nirdatafitreslist = [os.path.expandvars(filepath) for filepath in self.nirdatafitreslist]

		self.opticaldatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optical.FITRES.TEXT',
									  '$RAISIN_ROOT/cosmo/output/fit_optical/CfA_RAISIN_optical.FITRES.TEXT',
									  '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optical.FITRES.TEXT',
									  '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optical.FITRES.TEXT']
		self.opticaldatafitreslist = [os.path.expandvars(filepath) for filepath in self.opticaldatafitreslist]

		
		self.outfitres = 'output/fitres/RAISIN_stat.fitres'
		self.figdir = 'figs'
		
	def add_options(self):
		pass

	def mkcuts(self):
		pass
	
	def mk_biascor_files(self):
		pass

	def apply_biascor(self):
		# can interpolate for each SN individually because samples are so small
		pass

	def mk_sim_validplots(self):
		# make sure sim/data line up for all three sims
		# will have to worry about systematics down the road
		
		for simfitreslist,datafitreslist,label in zip([self.nirsimfitreslist,self.opticalsimfitreslist],
													  [self.nirdatafitreslist,self.opticaldatafitreslist],
													  ['NIR','Optical']):
			for i,survey in enumerate(['CSP','CfA','PS1','DES']):
				plt.clf()
				hist = ovdatamc.ovhist()
				parser = hist.add_options(usage='')
				options,  args = parser.parse_args()
				hist.options = options
				hist.options.journal = True
				hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.7,1.3)]
				hist.options.nbins = 10
				hist.options.clobber = True
				hist.options.outfile = 'figs/sim_%s_%s.png'%(survey,label)

				hist.options.histvar = ['MURES','SNRMAX1','AV','STRETCH']
				print(datafitreslist[i],simfitreslist[i])
				hist.main(datafitreslist[i],glob.glob('%s/*/FITOPT000.FITRES'%simfitreslist[i])[0])
				
	def mk_biascor_validplots(self):
		# make biascor plots
		pass

class sim:
	def __init__(self):
		pass

	def runsim(self):
		# LOW-Z
		os.system('sim_SNmix.pl $RAISIN_ROOT/cosmo/sim/inputs/lowz/SIMGEN_MASTER_LOWZ.INPUT')

		# PS1
		os.system('sim_SNmix.pl $RAISIN_ROOT/cosmo/sim/inputs/PS1/SIMGEN_MASTER_PS1SPEC.INPUT')

		# DES
		os.system('sim_SNmix.pl $RAISIN_ROOT/cosmo/sim/inputs/DES/SIMGEN_MASTER_DESSPEC.INPUT')
		
	def runfit_sim(self):
		# optical
		os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/CSP_RAISIN_optnir.nml')
		os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/CfA_RAISIN_optnir.nml')
		os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/PS1_RAISIN_optnir.nml')
		os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/DES_RAISIN_optnir.nml')
		
		# NIR
		os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/CSP_RAISIN_NIR.nml')
		os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/CfA_RAISIN_NIR.nml')
		os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/PS1_RAISIN_NIR.nml')
		os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/DES_RAISIN_NIR.nml')

class lcfit:
	def __init__(self):
		self.CSP_optical_fitresfile = 'output/fit_nir/CSP_RAISIN_optical.FITRES.TEXT'
		self.CfA_optical_fitresfile = 'output/fit_nir/CfA_RAISIN_optical.FITRES.TEXT'

		self.nirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN.FITRES.TEXT',
								  '$RAISIN_ROOT/cosmo/output/fit_nir/CfA_RAISIN.FITRES.TEXT',
								  '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN.FITRES.TEXT',
								  '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN.FITRES.TEXT']
		self.nirdatafitreslist = [os.path.expandvars(filepath) for filepath in self.nirdatafitreslist]

		self.opticaldatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optical.FITRES.TEXT',
									  '$RAISIN_ROOT/cosmo/output/fit_optical/CfA_RAISIN_optical.FITRES.TEXT',
									  '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optical.FITRES.TEXT',
									  '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optical.FITRES.TEXT']
		self.opticaldatafitreslist = [os.path.expandvars(filepath) for filepath in self.opticaldatafitreslist]

		
	def runfit_data(self):
		# optical
		os.system('snlc_fit.exe fit/CSP_RAISIN_optnir.nml')
		os.system('snlc_fit.exe fit/CfA_RAISIN_optnir.nml')
		os.system('snlc_fit.exe fit/PS1_RAISIN_optnir.nml')
		os.system('snlc_fit.exe fit/DES_RAISIN_optnir.nml')

		# nir
		os.system('snlc_fit.exe fit/CSP_RAISIN.nml')
		os.system('snlc_fit.exe fit/CfA_RAISIN.nml')
		os.system('snlc_fit.exe fit/PS1_RAISIN.nml')
		os.system('snlc_fit.exe fit/DES_RAISIN.nml')

	def plot_data(self):
		os.system('python raisin_cosmo/plot_snana.py -v CSPDR3 -o figs -a -20 -f fit/CSP_RAISIN_optnir.nml --plotAll -F output/fit_optical/CSP_RAISIN_optical.FITRES.TEXT')
		os.system('python raisin_cosmo/plot_snana.py -v CfAIR2 -o figs -a -20 -f fit/CfA_RAISIN_optnir.nml --plotAll -F output/fit_optical/CfA_RAISIN_optical.FITRES.TEXT')
		os.system('python raisin_cosmo/plot_snana.py -v PS1_RAISIN -o figs -a -20 -f fit/PS1_RAISIN_optnir.nml --plotAll -F output/fit_optical/CSP_RAISIN_optical.FITRES.TEXT')
		os.system('python raisin_cosmo/plot_snana.py -v DES_RAISIN -o figs -a -20 -f fit/DES_RAISIN_optnir.nml --plotAll -F output/fit_optical/CSP_RAISIN_optical.FITRES.TEXT')
		
	def add_pkmjd(self):
		CSPDR3_files = glob.glob('data/Photometry/CSPDR3/*.DAT')
		self.csp_opt_fr = txtobj(self.CSP_optical_fitresfile,fitresheader=True)

		for c in CSPDR3_files:
			if 'PKMJD' in c: continue
			with open(c) as fin:
				with open(c.replace('.DAT','.PKMJD.DAT'),'w') as fout:
					for line in fin:
						if line.startswith('PEAKMJD'):
							try:
								print('PEAKMJD:  %.2f                # from optical SNooPy fit'%(
									self.csp_opt_fr.PKMJD[self.csp_opt_fr.CID == SNID][0]),file=fout)
							except:
								sn = snana.SuperNova(c)
								print(SNID,sn.REDSHIFT_FINAL)
								print(line.replace('\n',''),file=fout)
						elif line.startswith('SNID'):
							SNID = line.split()[1].replace('\n','')
							print(line.replace('\n',''),file=fout)
						else:
							print(line.replace('\n',''),file=fout)

		# need to do the same for CfA


	def hubbleplot():
		pass
		
def errfnc(x):
    return(np.std(x)/np.sqrt(len(x)))

def main(frsim_ps1='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
		 frsim_ps1_oir='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY_OIR/FITOPT000.FITRES',
		 #frsim_ps1='output/sim/PS1_RAISIN_SIM_NIR.FITRES.TEXT',
		 #output/sim/PS1_RAISIN_SIM_NIR/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
		 frsim_des='output/sim/DES_RAISIN_SIM.FITRES.TEXT'): #'output/sim/DES_RAISIN_SIM_NIR.FITRES.TEXT'):

	plt.clf()
	ax = plt.axes() #subplot(121)
	
	frp = txtobj(frsim_ps1,fitresheader=True)
	frpo = txtobj(frsim_ps1_oir,fitresheader=True)
	frd = txtobj(frsim_des,fitresheader=True)

	zbins = np.linspace(0.01,0.7,20)
	#iSNR = frp.SNRMAX1 > 50
	#ps1_bin = binned_statistic(
	#	frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
	#	statistic='median',bins=zbins).statistic
	#ps1_errbin = binned_statistic(
	#	frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
	#	statistic=errfnc,bins=zbins).statistic
	#ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
	#			fmt='o-',color='b',capsize=0,lw=2)
	ps1_bin = binned_statistic(
		frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic='median',bins=zbins).statistic
	ps1_bin_oir = binned_statistic(
		frpo.zCMB,(frpo.mB-frpo.SIM_mB)+0.14*(frpo.x1-frpo.SIM_x1) - 3.1*(frpo.c-frpo.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic='median',bins=zbins).statistic

	#ps1_bin = binned_statistic(
	#	frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
	#	statistic='median',bins=zbins).statistic
	ps1_errbin = binned_statistic(
		frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic=errfnc,bins=zbins).statistic
	ps1_errbin_oir = binned_statistic(
		frpo.zCMB,(frpo.mB-frpo.SIM_mB)+0.14*(frpo.x1-frpo.SIM_x1) - 3.1*(frpo.c-frpo.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic=errfnc,bins=zbins).statistic

	#ps1_errbin = binned_statistic(
	#	frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
	#	statistic=errfnc,bins=zbins).statistic

	
	#des_bin = binned_statistic(
	#	frd.zCMB,frd.DLMAG-frd.SIM_DLMAG,
	#	statistic='median',bins=zbins).statistic

	#ax2 = plt.subplot(122)
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
				fmt='o-',color='r',capsize=0,lw=2)
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin_oir,yerr=ps1_errbin_oir,
				fmt='o-',color='b',capsize=0,lw=2)

	#ax.errorbar((zbins[1:]+zbins[:-1])/2.,des_bin,
	#			fmt='o-',color='b',capsize=0,lw=2)
	ax.set_xlabel('$z_{CMB}$',fontsize=15)
	ax.set_ylabel('$\mu - \mu_{sim}$',fontsize=15)
	import pdb; pdb.set_trace()
	
	return

def smearmodelcomp(
		frsim_ps1='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
		frsim_ps1_oir='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY_OIR_NEWSIMLIB_z/FITOPT000.FITRES',
		#frsim_ps1='output/sim/PS1_RAISIN_SIM_NIR.FITRES.TEXT',
		#output/sim/PS1_RAISIN_SIM_NIR/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
		frsim_des='output/sim/DES_RAISIN_SIM.FITRES.TEXT'): #'output/sim/DES_RAISIN_SIM_NIR.FITRES.TEXT'):

	plt.clf()
	ax = plt.axes() #subplot(121)
	
	frp = txtobj(frsim_ps1,fitresheader=True)
	frpo = txtobj(frsim_ps1_oir,fitresheader=True)
	frd = txtobj(frsim_des,fitresheader=True)

	zbins = np.linspace(0.01,0.7,20)
	#iSNR = frp.SNRMAX1 > 50
	#ps1_bin = binned_statistic(
	#	frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
	#	statistic='median',bins=zbins).statistic
	#ps1_errbin = binned_statistic(
	#	frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
	#	statistic=errfnc,bins=zbins).statistic
	#ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
	#			fmt='o-',color='b',capsize=0,lw=2)
	ps1_bin = binned_statistic(
		frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic='median',bins=zbins).statistic
	ps1_bin_oir = binned_statistic(
		frpo.zCMB,(frpo.mB-frpo.SIM_mB)+0.14*(frpo.x1-frpo.SIM_x1) - 3.1*(frpo.c-frpo.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic='median',bins=zbins).statistic

	#ps1_bin = binned_statistic(
	#	frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
	#	statistic='median',bins=zbins).statistic
	ps1_errbin = binned_statistic(
		frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic=errfnc,bins=zbins).statistic
	ps1_errbin_oir = binned_statistic(
		frpo.zCMB,(frpo.mB-frpo.SIM_mB)+0.14*(frpo.x1-frpo.SIM_x1) - 3.1*(frpo.c-frpo.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic=errfnc,bins=zbins).statistic
	
	
	#ps1_errbin = binned_statistic(
	#	frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
	#	statistic=errfnc,bins=zbins).statistic

	
	#des_bin = binned_statistic(
	#	frd.zCMB,frd.DLMAG-frd.SIM_DLMAG,
	#	statistic='median',bins=zbins).statistic

	#ax2 = plt.subplot(122)
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
				fmt='o-',color='r',capsize=0,lw=2,label='G10 model')
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin_oir,yerr=ps1_errbin_oir,
				fmt='o-',color='b',capsize=0,lw=2,label='SNooPy OIR model')
	ax.legend()
	#ax.errorbar((zbins[1:]+zbins[:-1])/2.,des_bin,
	#			fmt='o-',color='b',capsize=0,lw=2)
	ax.set_xlabel('$z_{CMB}$',fontsize=15)
	ax.set_ylabel('$\mu - \mu_{sim}$',fontsize=15)
	import pdb; pdb.set_trace()
	
	return

def saltvsnoopy(
		frsim_ps1='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY_OIR_newsimlib_z/FITOPT000.FITRES',
		frsim_ps1_oir='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY_OIR_real_nironly/FITOPT000.FITRES',
		#frsim_ps1='output/sim/PS1_RAISIN_SIM_NIR.FITRES.TEXT',
		#output/sim/PS1_RAISIN_SIM_NIR/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
		frsim_des='output/sim/DES_RAISIN_SIM.FITRES.TEXT'): #'output/sim/DES_RAISIN_SIM_NIR.FITRES.TEXT'):

	plt.clf()
	ax = plt.axes() #subplot(121)
	
	frp = txtobj(frsim_ps1,fitresheader=True)
	frpo = txtobj(frsim_ps1_oir,fitresheader=True)
	frd = txtobj(frsim_des,fitresheader=True)
	#igood = frpo.AVERR < 0.05
	#for k in frpo.__dict__.keys():
	#	frpo.__dict__[k] = frpo.__dict__[k][igood]
	
	zbins = np.linspace(0.01,0.7,20)
	#iSNR = frp.SNRMAX1 > 50
	#ps1_bin = binned_statistic(
	#	frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
	#	statistic='median',bins=zbins).statistic
	#ps1_errbin = binned_statistic(
	#	frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
	#	statistic=errfnc,bins=zbins).statistic
	#ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
	#			fmt='o-',color='b',capsize=0,lw=2)
	ps1_bin = binned_statistic(
		frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic='median',bins=zbins).statistic
	ps1_bin_oir = binned_statistic(
		frpo.zCMB,frpo.DLMAG-frpo.SIM_DLMAG,
		statistic='median',bins=zbins).statistic
	ps1_stretch_oir = binned_statistic(
		frpo.zCMB,frpo.STRETCH-frpo.SIM_STRETCH,
		statistic='median',bins=zbins).statistic
	ps1_av_oir = binned_statistic(
		frpo.zCMB,frpo.AV-frpo.SIM_AV,
		statistic='median',bins=zbins).statistic

	
	#ps1_bin = binned_statistic(
	#	frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
	#	statistic='median',bins=zbins).statistic
	ps1_errbin = binned_statistic(
		frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic=errfnc,bins=zbins).statistic
	ps1_errbin_oir = binned_statistic(frpo.zCMB,
		frpo.DLMAG - frpo.SIM_DLMAG,
		statistic=errfnc,bins=zbins).statistic
	
	
	#ps1_errbin = binned_statistic(
	#	frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
	#	statistic=errfnc,bins=zbins).statistic

	
	#des_bin = binned_statistic(
	#	frd.zCMB,frd.DLMAG-frd.SIM_DLMAG,
	#	statistic='median',bins=zbins).statistic

	#ax2 = plt.subplot(122)
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
				fmt='o-',color='r',capsize=0,lw=2,label='G10 model')
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin_oir,yerr=ps1_errbin_oir,
				fmt='o-',color='b',capsize=0,lw=2,label='SNooPy OIR model')
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_stretch_oir,yerr=ps1_errbin_oir,
				fmt='o-',color='C0',capsize=0,lw=2,label='SNooPy OIR model - $s$')
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_av_oir,yerr=ps1_errbin_oir,
				fmt='o-',color='C1',capsize=0,lw=2,label='SNooPy OIR model - A_V')

	ax.legend()
	#ax.errorbar((zbins[1:]+zbins[:-1])/2.,des_bin,
	#			fmt='o-',color='b',capsize=0,lw=2)
	ax.set_xlabel('$z_{CMB}$',fontsize=15)
	ax.set_ylabel('$\mu - \mu_{sim}$',fontsize=15)
	import pdb; pdb.set_trace()
	
	return

def saltvsnoopy_des(
		frsim_des='output/sim/DES_RAISIN_SIM/DES_RAISIN_SALT2/FITOPT000.FITRES',
		frsim_des_oir='output/sim/DES_RAISIN_SIM_NIR/DES_RAISIN_SNOOPY_OIR/FITOPT000.FITRES'):
	plt.clf()
	ax = plt.axes() #subplot(121)
	
	frp = txtobj(frsim_des,fitresheader=True)
	frpo = txtobj(frsim_des_oir,fitresheader=True)
	frd = txtobj(frsim_des,fitresheader=True)
	#igood = frpo.AVERR < 0.05
	#for k in frpo.__dict__.keys():
	#	frpo.__dict__[k] = frpo.__dict__[k][igood]
	
	zbins = np.linspace(0.01,0.7,20)
	#iSNR = frp.SNRMAX1 > 50
	#ps1_bin = binned_statistic(
	#	frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
	#	statistic='median',bins=zbins).statistic
	#ps1_errbin = binned_statistic(
	#	frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
	#	statistic=errfnc,bins=zbins).statistic
	#ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
	#			fmt='o-',color='b',capsize=0,lw=2)
	ps1_bin = binned_statistic(
		frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic='median',bins=zbins).statistic
	ps1_bin_oir = binned_statistic(
		frpo.zCMB,frpo.DLMAG-frpo.SIM_DLMAG,
		statistic='median',bins=zbins).statistic
	ps1_stretch_oir = binned_statistic(
		frpo.zCMB,frpo.STRETCH-frpo.SIM_STRETCH,
		statistic='median',bins=zbins).statistic
	ps1_av_oir = binned_statistic(
		frpo.zCMB,frpo.AV-frpo.SIM_AV,
		statistic='median',bins=zbins).statistic

	
	#ps1_bin = binned_statistic(
	#	frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
	#	statistic='median',bins=zbins).statistic
	ps1_errbin = binned_statistic(
		frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
		statistic=errfnc,bins=zbins).statistic
	ps1_errbin_oir = binned_statistic(frpo.zCMB,
		frpo.DLMAG - frpo.SIM_DLMAG,
		statistic=errfnc,bins=zbins).statistic
	
	
	#ps1_errbin = binned_statistic(
	#	frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
	#	statistic=errfnc,bins=zbins).statistic

	
	#des_bin = binned_statistic(
	#	frd.zCMB,frd.DLMAG-frd.SIM_DLMAG,
	#	statistic='median',bins=zbins).statistic

	#ax2 = plt.subplot(122)
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
				fmt='o-',color='r',capsize=0,lw=2,label='G10 model')
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin_oir,yerr=ps1_errbin_oir,
				fmt='o-',color='b',capsize=0,lw=2,label='SNooPy OIR model')
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_stretch_oir,yerr=ps1_errbin_oir,
				fmt='o-',color='C0',capsize=0,lw=2,label='SNooPy OIR model - $s$')
	ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_av_oir,yerr=ps1_errbin_oir,
				fmt='o-',color='C1',capsize=0,lw=2,label='SNooPy OIR model - A_V')

	ax.legend()
	#ax.errorbar((zbins[1:]+zbins[:-1])/2.,des_bin,
	#			fmt='o-',color='b',capsize=0,lw=2)
	ax.set_xlabel('$z_{CMB}$',fontsize=15)
	ax.set_ylabel('$\mu - \mu_{sim}$',fontsize=15)
	import pdb; pdb.set_trace()
	
	return

if __name__ == "__main__":
	#main()
	#smearmodelcomp()
	#saltvsnoopy()
	#saltvsnoopy_des()

	sm = sim()
	sm.runsim()
	#sm.runfit()
	
	lcf = lcfit()
	lcf.add_pkmjd()
	
	#bc = biascor()
	#bc.mk_sim_validplots()
	
	#bc.mkcuts()
	#bc.mk_biascor_files()
	#bc.apply_biascor()
