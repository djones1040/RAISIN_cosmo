import pylab as plt
plt.ion()
import numpy as np
import cosmo
from txtobj import txtobj
from astropy.stats import sigma_clipped_stats
import copy
#https://jiffyclub.github.io/palettable/

fitresfile_optical_ps1 = 'output/fit_optical/PS1_RAISIN_optical.FITRES.TEXT'
fitresfile_optical_des = 'output/fit_optical/DES_RAISIN_optical.FITRES.TEXT'
fitresfile_optical_lowz = 'output/fit_optical/LOWZ_RAISIN_optical.FITRES.TEXT'
fitresfile_ps1 = 'output/fit_nir/PS1_RAISIN.FITRES.TEXT'
fitresfile_des = 'output/fit_nir/DES_RAISIN.FITRES.TEXT'
fitresfile_lowz = 'output/fit_nir/LOWZ_RAISIN.FITRES.TEXT'

distkeys = ['Ymu','Jmu','Hmu']
disterrkeys = ['Ymuerr','Jmuerr','Hmuerr']

GOOD_RAISIN_CIDS = ['PScH540118', 'PScH540087', 'DES16C3cmy', 'DES15E2nlz',
					'PScG530251', 'DES16S2afz', 'PScJ550202', 'PScJ440005',
					'PScK450339', 'PScA470110', 'DES16E2clk', #'PScF510457',
					'DES16E1dcx', 'DES15X2nkz', 'DES15X2kvt', 'DES16S1agd',
					'DES16E2cqq', 'PScB480464', 'DES15C1nhv',
					'DES16S1bno',
					'SNABELL370', 'PScD500100', 'PScD500301', 'PScF520062',
					'PScC490037', 'DES16C2cva', 'PScJ560027', #'PScF520107',
					'DES16C1cim', 'PScA470240', 'PScF520188', 'DES16X3zd',
					'DES15X2mey', 'PScB480794', 'DES15E2mhy', 'DES16X3cry',
					'PScJ560054', 'PScJ440236', 'DES15C3odz', 'PScC490521',
					'DES15E2uc']

def add_optical_info(fr):
	opt_fr = txtobj(fitresfile_optical_ps1,fitresheader=True)
	opt_fr2 = txtobj(fitresfile_optical_des,fitresheader=True)
	opt_fr3 = txtobj(fitresfile_optical_lowz,fitresheader=True)
	for k in opt_fr.__dict__.keys():
		if k != 'MWEBV':
			opt_fr.__dict__[k] = np.concatenate((
				opt_fr.__dict__[k],opt_fr2.__dict__[k],opt_fr3.__dict__[k],))
	
	fr.opt_stretch,fr.opt_ebv,fr.opt_fitprob = \
		np.array([-99.0]*len(fr.CID)),np.array([-99.0]*len(fr.CID)),np.array([-99.0]*len(fr.CID))

	for j,i in enumerate(fr.CID):
		if i == 'SNABELL370': fr.opt_stretch[j] = 1; fr.opt_ebv[j] = 0
		else:
			if not i in opt_fr.CID:
				print(i)#raise RuntimeError('missing optical fit for %s'%fr.CID[j])
				continue
			fr.opt_stretch[j] = opt_fr.STRETCH[opt_fr.CID == i][0]
			fr.opt_ebv[j] = opt_fr.AV[opt_fr.CID == i][0]/3.1
			fr.opt_fitprob[j] = opt_fr.FITPROB[opt_fr.CID == i][0]

	iBad = (fr.opt_ebv > 0.3) | (fr.opt_stretch < 0.7) | (fr.opt_stretch > 1.3) #| (fr.opt_fitprob < 0.001)
	print(fr.CID[iBad],fr.opt_stretch[iBad],fr.opt_ebv[iBad])

	iGood = (fr.opt_ebv < 0.3) & (fr.opt_stretch > 0.7) & (fr.opt_stretch < 1.3) #& (fr.opt_fitprob > 0.001)
	for k in fr.__dict__.keys():
		fr.__dict__[k] = fr.__dict__[k][iGood]

	return fr
		
def hubbleplot():
	plt.rcParams['figure.figsize'] = (8,4)
	plt.rcParams['font.size'] = 13

	plt.clf()
	plt.subplots_adjust(left=None, bottom=None, right=None,
						top=None, wspace=0, hspace=0)
	
	ax1p1 = plt.axes([0.1,0.45,0.33,0.45])
	ax1p2 = plt.axes([0.48,0.45,0.33,0.45],sharey=ax1p1)
	ax2p1 = plt.axes([0.1,0.13,0.33,0.32])
	ax2p2 = plt.axes([0.48,0.13,0.33,0.32],sharey=ax2p1)
	ax3 = plt.axes([0.81,0.13,0.1,0.32],sharey=ax2p1)
	ax3.yaxis.tick_right()
	ax3.set_xticks([])

	import palettable

	from palettable.colorbrewer.qualitative import Dark2_8 as palettable_color
	for ax in [ax1p1,ax1p2,ax2p1,ax2p2,ax3]:
		ax.set_prop_cycle('color', palettable_color.mpl_colors)

	ax1p1.spines['right'].set_visible(False)
	ax1p2.spines['left'].set_visible(False)
	ax2p1.spines['right'].set_visible(False)
	ax2p2.spines['left'].set_visible(False)

	# get the distance moduli
	fr = txtobj(fitresfile_ps1,fitresheader=True)
	fr.SURVEY = np.array(['RAISIN1']*len(fr.CID))
	fr2 = txtobj(fitresfile_des,fitresheader=True)
	fr2.SURVEY = np.array(['RAISIN2']*len(fr2.CID))
	for k in fr.__dict__.keys():
		fr.__dict__[k] = np.append(fr.__dict__[k],fr2.__dict__[k])

	print('%i RAISIN SNe before cuts'%len(fr.zCMB))
	fr = mkcuts(fr,raisincuts=True)
	print('%i RAISIN SNe after cuts'%len(fr.zCMB))
	fr = add_optical_info(fr)
	print('%i RAISIN SNe after optical cuts'%len(fr.zCMB))

	frlowz = txtobj(fitresfile_lowz,fitresheader=True)
	print('%i low-z SNe before cuts'%len(frlowz.zCMB))
	frlowz = mkcuts(frlowz)
	print('%i low-z SNe after cuts'%len(frlowz.zCMB))
	frlowz = add_optical_info(frlowz)
	print('%i low-z SNe after optical cuts'%len(frlowz.zCMB))

	avgresid = sigma_clipped_stats(
		np.append(fr.DLMAG-cosmo.mu(fr.zCMB),
				  frlowz.DLMAG-cosmo.mu(frlowz.zCMB)))[1]
	fr.muresid = fr.DLMAG - cosmo.mu(fr.zCMB) - avgresid
	frlowz.muresid = frlowz.DLMAG - cosmo.mu(frlowz.zCMB) - avgresid
	fr.DLMAG -= avgresid
	#fr.muavgresid -= avgresid
	frlowz.DLMAG -= avgresid
	#frlowz.muavgresid -= avgresid

	# plot the hubble resids
	iR1 = fr.SURVEY == 'RAISIN1'
	iR2 = fr.SURVEY == 'RAISIN2'
	ax1p1.errorbar(frlowz.zCMB,frlowz.DLMAG,yerr=frlowz.DLMAGERR,fmt='o',
				   label='low-$z$',color=palettable_color.mpl_colors[0])
	ax1p1.errorbar(fr.zCMB[iR1],fr.DLMAG[iR1],yerr=fr.DLMAGERR[iR1],fmt='o',
				   label='RAISIN1',color=palettable_color.mpl_colors[1])
	ax1p1.errorbar(fr.zCMB[iR2],fr.DLMAG[iR2],yerr=fr.DLMAGERR[iR2],fmt='o',
				   label='RAISIN2',color=palettable_color.mpl_colors[2])

	ax1p2.errorbar(fr.zCMB[iR1],fr.DLMAG[iR1],yerr=fr.DLMAGERR[iR1],fmt='o',
				   label='RAISIN1',color=palettable_color.mpl_colors[1])
	ax1p2.errorbar(fr.zCMB[iR2],fr.DLMAG[iR2],yerr=fr.DLMAGERR[iR2],fmt='o',
				   label='RAISIN2',color=palettable_color.mpl_colors[2])
	plotlcdm(ax1p1)
	plotlcdm(ax1p2)
	ax2p1.errorbar(frlowz.zCMB,frlowz.muresid,yerr=frlowz.DLMAGERR,fmt='o',
				   label='low-$z$',color=palettable_color.mpl_colors[0])
	ax2p2.errorbar(fr.zCMB[iR1],fr.muresid[iR1],yerr=fr.DLMAGERR[iR1],fmt='o',
				   color=palettable_color.mpl_colors[1])
	ax2p2.errorbar(fr.zCMB[iR2],fr.muresid[iR2],yerr=fr.DLMAGERR[iR2],fmt='o',
				   color=palettable_color.mpl_colors[2])
	ax2p1.axhline(0,color='k',lw=2)
	ax2p2.axhline(0,color='k',lw=2)
	print('N_RAISIN = %i'%len(fr.zCMB))
	print('N_LOWZ = %i'%len(frlowz.zCMB))
	# making the plots look good
	ax1p1.set_xlim([0.01,0.05])
	ax1p2.set_xlim([0.2,0.7])
	ax2p1.set_xlim([0.01,0.05])
	ax2p2.set_xlim([0.2,0.7])
	ax1p1.legend(loc='upper left')

	ax3.hist(np.append(fr.muresid,frlowz.muresid[frlowz.zCMB > 0.01]),bins=np.arange(-1,1,0.1),
			 orientation='horizontal',color='r',density=True,alpha=0.7,ec='r',lw=2,label='NIR')
	frpan = txtobj('hlsp_ps1cosmo_panstarrs_gpc1_all_model_v1_ancillary-g10.fitres.txt',fitresheader=True)
	ax3.hist(frpan.MURES,bins=np.arange(-1,1,0.1),
			 orientation='horizontal',color='b',density=True,
			 alpha=1.0,ec='k',lw=2,label='Optical\n(Pantheon)',histtype='step')
	ax3.legend(prop={'size':10},bbox_to_anchor=(0.0,1.6),bbox_transform=ax3.transAxes,loc='upper left')

	from matplotlib.ticker import NullFormatter
	for ax in [ax1p1,ax1p2]:
		ax.set_xscale('log')
		ax.xaxis.set_major_formatter(NullFormatter())
		ax.xaxis.set_minor_formatter(NullFormatter())

		ax.set_xticks([])
		ax.set_xlabel(r'$z_{\mathrm{CMB}}$',fontsize=15,labelpad=0)

		ax.set_ylabel('$\mu$ (mag)',fontsize=15)
		ax.set_ylim([33,45])
		ax.tick_params(top="off",bottom="on",direction="inout",length=8, width=2)
	for ax in [ax2p1,ax2p2]:
		ax.set_xscale('log')
		ax.xaxis.set_major_formatter(NullFormatter())
		ax.xaxis.set_minor_formatter(NullFormatter())

		ax.set_xlabel(r'$z_{\mathrm{CMB}}$',fontsize=15,labelpad=0)

		ax.set_ylabel('$\mu - \mu_{\Lambda CDM}$ (mag)',fontsize=15,labelpad=0)
		ax.set_ylim([-1,1])
		ax.set_ylim([-0.5,0.5])
	
	ax.tick_params(top="off",bottom="on",direction="inout",length=8, width=2)

	ax1p2.set_ylabel('')
	ax2p2.set_ylabel('')
	ax1p1.set_xticks([0.01,0.03])
	ax1p1.set_xticklabels(['0.01','0.03'])
	ax1p2.set_xticks([0.3,0.5,0.7])
	ax1p2.set_xticklabels(['0.3','0.5','0.7'])
	ax2p1.set_xticks([0.01,0.03])
	ax2p1.set_xticklabels(['0.01','0.03'])
	ax2p2.set_xticks([0.3,0.5,0.7])
	ax2p2.set_xticklabels(['0.3','0.5','0.7'])
	ax1p2.set_xlim([0.2,0.8])
	ax2p2.set_xlim([0.2,0.8])
	ax1p1.set_yticks([35,37.5,40,42.5,45])
	ax1p1.set_xticklabels([])
	
	ax1p2.yaxis.tick_right()
	plt.setp(ax1p2.get_yticklabels(), visible=False)
	
	ax2p2.yaxis.tick_right()
	plt.setp(ax2p2.get_yticklabels(), visible=False)


	
	d = .015 # how big to make the diagonal lines in axes coordinates
	# arguments to pass plot, just so we don't keep repeating them
	kwargs = dict(transform=ax1p1.transAxes, color='k', clip_on=False)
	ax1p1.plot((1-d,1+d), (-d,+d), **kwargs)
	ax1p1.plot((1-d,1+d),(1-d,1+d), **kwargs)

	kwargs.update(transform=ax1p2.transAxes)  # switch to the bottom axes
	ax1p2.plot((-d,+d), (1-d,1+d), **kwargs)
	ax1p2.plot((-d,+d), (-d,+d), **kwargs)

	kwargs = dict(transform=ax2p1.transAxes, color='k', clip_on=False)
	ax2p1.plot((1-d,1+d), (-d,+d), **kwargs)
	ax2p1.plot((1-d,1+d),(1-d,1+d), **kwargs)

	kwargs.update(transform=ax2p2.transAxes)  # switch to the bottom axes
	ax2p2.plot((-d,+d), (1-d,1+d), **kwargs)
	ax2p2.plot((-d,+d), (-d,+d), **kwargs)

	print(fr.CID)
	
	import pdb; pdb.set_trace()
	
def mkcuts(fr,raisincuts=False):

	if 'ebv' in fr.__dict__.keys():
		iCut = fr.ebv < 0.25
		for k in fr.__dict__.keys():
			fr.__dict__[k] = fr.__dict__[k][iCut]

	if raisincuts:
		good = np.array([],dtype=int)
		for j,i in enumerate(fr.CID):
			if i in GOOD_RAISIN_CIDS:
				good = np.append(good,j)
		for k in fr.__dict__.keys():
			fr.__dict__[k] = fr.__dict__[k][good]
			
	iCut = (fr.CID != 'snf20080514-002') & (fr.CID != 'sn2005bo')
	for k in fr.__dict__.keys():
		fr.__dict__[k] = fr.__dict__[k][iCut]

	iCut = fr.DLMAGERR < 0.4
	for k in fr.__dict__.keys():
		fr.__dict__[k] = fr.__dict__[k][iCut]

	return fr
	
def hubbleplotpars(ax):

	ax.set_xscale('log')
	ax.set_xlabel(r'$z_{\mathrm{CMB}}$',fontsize=15,labelpad=0)
	ax.set_xlim([0.01,0.75])

	ax.set_ylabel('$\mu$ (mag)',fontsize=15)
	ax.set_ylim([33,45])
	ax.tick_params(top="off",bottom="on",direction="inout",length=8, width=2)

	ax.set_xticklabels([])

def hubbleresidplotpars(ax):

	ax.set_xscale('log')
	ax.set_xlabel(r'$z_{\mathrm{CMB}}$',fontsize=15,labelpad=0)
	ax.set_xlim([0.01,0.75])

	ax.set_ylabel('$\mu - \mu_{\Lambda CDM}$ (mag)',fontsize=15,labelpad=0)
	ax.set_ylim([-1,1])
	ax.set_ylim([-0.5,0.5])

	ax.tick_params(top="off",bottom="on",direction="inout",length=8, width=2)

	
def plotlcdm(ax):

	zrange = np.arange(0,1,0.01)
	ax.plot(zrange,cosmo.mu(zrange),color='k',lw=2)

def weighted_avg_and_err(values, weights):
	"""
	Return the weighted average and standard deviation.
	values, weights -- Numpy ndarrays with the same shape.
	"""
	average = np.average(values, weights=weights)
	variance = np.average((values-average)**2, weights=weights)	 # Fast and numerically precise
	return (average, np.sqrt(variance/len(values)))			

if __name__ == "__main__":
	hubbleplot()
