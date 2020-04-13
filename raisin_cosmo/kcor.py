#!/usr/bin/env python
# D. Jones - 11/21/19
"""see how different SN Ia models affect distance measurements via k-corrections"""

import os
import numpy as np
from txtobj import txtobj
from util.common import *
import pylab as plt
plt.ion()

class kcor:

	def __init__(self):

		self.snsed = ['kcor/snsed/Hsiao07.dat',
					  'kcor/snsed/sn91t_flux.v1.1.dat',
					  #'kcor/snsed/salt2_m0_P18.dat',
					  'kcor/snsed/snflux_1a_Nugent2002.dat']

		self.kcor_input_ps1 = 'kcor/kcor_PS1MD_NIR.input'
		self.kcor_input_des = 'kcor/kcor_DES_NIR.input'
		self.kcor_input_hst = ''

		self.nml_ps1 = 'fit/PS1_RAISIN.nml'
		self.nml_des = 'fit/DES_RAISIN.nml'
		self.nml_hst = ''
		
	def mk_kcor(self):
		for ia_template in self.snsed[1:2]:
			ia_template_name = ia_template.split('/')[-1].split('.')[0]

			for kcor in [self.kcor_input_ps1,self.kcor_input_des]:
				with open(kcor,'r') as fin:
					lines = fin.readlines()
				new_kcor = kcor.replace('.input','_%s.input'%ia_template_name)
				with open(new_kcor,'w') as fout:
					for line in lines:
						if not line.startswith('SN_SED') and not line.startswith('OUTFILE'):
							print(line.replace('\n',''),file=fout)
						elif line.startswith('SN_SED'):
							print('SN_SED: %s'%ia_template,file=fout)
						elif line.startswith('OUTFILE'):
							print('OUTFILE: %s'%kcor.replace('.input','_%s.fits'%ia_template_name),file=fout)

				os.system('kcor.exe %s'%new_kcor)
			
	def run_lcfit(self):
		for ia_template in self.snsed:
			ia_template_name = ia_template.split('/')[-1].split('.')[0]
			for nml in [self.nml_ps1,self.nml_des]:
				with open(nml,'r') as fin:
					lines = fin.readlines()
				new_nml = nml.replace('.nml','_%s.nml'%ia_template_name)
				with open(new_nml,'w') as fout:
					for line in lines:
						if 'TEXTFILE_PREFIX' not in line and 'KCOR_FILE' not in line:
							print(line.replace('\n',''),file=fout)
						elif 'TEXTFILE_PREFIX' in line:
							print("	 TEXTFILE_PREFIX = 'output/kcor_sys/%s'"%nml.split('/')[-1].replace('.nml','_%s'%ia_template_name),file=fout)
						elif 'KCOR_FILE' in line:
							print("	 KCOR_FILE = %s"%line.split()[-1].replace('.fits','_%s.fits'%ia_template_name),file=fout)

				print('snlc_fit.exe %s'%new_nml)
				import pdb; pdb.set_trace()
				os.system('snlc_fit.exe %s'%new_nml)

	def comp_distances(self):
		hsiao_frfile_ps1 = 'output/kcor_sys/PS1_RAISIN_Hsiao07.FITRES.TEXT'; hsiao_frfile_des = 'output/kcor_sys/DES_RAISIN_Hsiao07.FITRES.TEXT'
		nugent_frfile_ps1 = 'output/kcor_sys/PS1_RAISIN_snflux_1a_Nugent2002.FITRES.TEXT'; nugent_frfile_des = 'output/kcor_sys/DES_RAISIN_snflux_1a_Nugent2002.FITRES.TEXT'
		p18_frfile_ps1 = 'output/kcor_sys/PS1_RAISIN_sn91t_flux.FITRES.TEXT'; p18_frfile_des = 'output/kcor_sys/DES_RAISIN_sn91t_flux.FITRES.TEXT'

		frbase = txtobj(hsiao_frfile_ps1,fitresheader=True); frbase2 = txtobj(hsiao_frfile_des,fitresheader=True)
		frn = txtobj(nugent_frfile_ps1,fitresheader=True); frn2 = txtobj(nugent_frfile_des,fitresheader=True)
		frp = txtobj(p18_frfile_ps1,fitresheader=True); frp2 = txtobj(p18_frfile_des,fitresheader=True)
		frbase = concat_simple(frbase,frbase2)
		frn = concat_simple(frn,frn2)
		frp = concat_simple(frp,frp2)
		
		frbase = mkraisincuts(frbase)
		frn = mkraisincuts(frn)
		frp = mkraisincuts(frp)
		
		iSort = np.argsort(frbase.zHD)
		frbase = sort_fitres(frbase,iSort); frn = sort_fitres(frn,iSort); frp = sort_fitres(frp,iSort)

		plt.clf()
		plt.plot(frbase.zCMB,frn.DLMAG-frbase.DLMAG,'o-',label='Nugent02 - Hsiao+07')
		#plt.plot(frbase.zCMB,frn.DLMAG-frbase.DLMAG,'o-',label='Nugent91T - Hsiao+07')
		plt.xlabel('$z_{CMB}$',fontsize=15)
		plt.ylabel('$\Delta\mu$',fontsize=15)
		plt.axhline(0,color='k',lw=2)
		plt.legend()
		import pdb; pdb.set_trace()
		
if __name__ == "__main__":

	kc = kcor()
	#kc.mk_kcor()
	#kc.run_lcfit()
	kc.comp_distances()
