#!/usr/bin/env python
# D. Jones - 2/17/21
# evaluate prospects for Roman SN survey distances

# goal is a plot of Hubble residual scatter as a
# function of (1) minimum phase and (2) cadence
# with YJ and YJH light curves.

# for this we want reasonable, small-ish simulations of Roman
# survey at z < 0.7 using SNooPy
# cadence: 3,5,7,9,11 days
# min phase, controlling for number of data points?
# Roman won't choose this anyways so who cares,
# but perhaps useful for future NIR data.

# (1) SNANA sim-input files, starting from nominal stuff that Justin was using
# (2) need different simlibs to set the cadence, so go into nominal simlib and
#     just edit the MJDs.  I don't think I care about taking root(N) into account
#     here.  Or I can just do a 1-day cadence and edit the LC files themselves.
# (3) set redshift to z < 0.7, just do YJH
# (4) fit the data, keep the FITRES files
# (5) make the plots

_lcdirs_cadence = []
_cadences = []
_lcdirs_minphase = []
_minphase = []
import os
import glob
import numpy as np
import pylab as plt
plt.ion()
from txtobj import txtobj
import cosmo
import snana
from astropy.stats import sigma_clipped_stats

def edit_simlib():

    count = 1
    with open('3tier_Survey.SIMLIB') as fin, open('3tier_Survey_oneday.SIMLIB','w') as fout:
        for line in fin:
            if line.startswith('S:'):
                mjd = float(line.split()[1])
                lineparts = line.split()
                for i in range(5):
                    lineparts[1] = str(mjd + i)
                    lineparts[2] = str(count)
                    count += 1
                    print(' '.join(lineparts),file=fout)
            else:
                print(line.replace('\n',''),file=fout)
    return

def edit_lcs():
    # function for making the lightcurve versions
    for cadence in [1,3,5,7,9,11,13]:
        if not os.path.exists(f'ROMAN_NIR_cadence{cadence:.0f}'):
            os.system(f'mkdir ROMAN_NIR_cadence{cadence:.0f}')
        lcfiles = glob.glob('/usr/local/SNDATA_ROOT/SCRATCH/ROMAN_NIR/ROMAN_NIR_SN*.DAT')

        for l in lcfiles:
            sn = snana.SuperNova(l.replace('/usr/local/SNDATA_ROOT/SCRATCH/ROMAN_NIR',
                                f'ROMAN_NIR_cadence{cadence:.0f}'))
            sn.FLT = sn.BAND
            nobs = len(sn.FLT)
            
            with open(l) as fin, \
                 open(l.replace('/usr/local/SNDATA_ROOT/SCRATCH/ROMAN_NIR',
                                f'ROMAN_NIR_cadence{cadence:.0f}'),'w') as fout:
                lastmjdZ,lastmjdY,lastmjdJ,lastmjdH = None,None,None,None
                for line in fin:
                    if line.startswith('NOBS:'):
                        print(f'NOBS: {nobs:.0f}',file=fout)
                    elif not line.startswith('OBS:'):
                        print(line.replace('\n',''),file=fout)
                    else:
                        mjd = float(line.split()[1])
                        filt = line.split()[2]
                        if filt == 'Z': lastmjd = lastmjdZ
                        elif filt == 'Y': lastmjd = lastmjdY
                        elif filt == 'J': lastmjd = lastmjdJ
                        elif filt == 'H': lastmjd = lastmjdH

                        if filt == 'Z' and lastmjd is None:
                            lastmjdZ = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'Z' and mjd >= lastmjdZ+cadence:
                            lastmjdZ = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'Y' and lastmjd is None:
                            lastmjdY = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'Y' and mjd >= lastmjdY+cadence:
                            lastmjdY = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'J' and lastmjd is None:
                            lastmjdJ = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'J' and mjd >= lastmjdJ+cadence:
                            lastmjdJ = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'H' and lastmjd is None:
                            lastmjdH = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'H' and mjd >= lastmjdH+cadence:
                            lastmjdH = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        
        os.chdir(f'ROMAN_NIR_cadence{cadence:.0f}')
        os.system(f'ls * > ROMAN_NIR_cadence{cadence:.0f}.LIST')
        with open(f'ROMAN_NIR_cadence{cadence:.0f}.README','w') as fout:
            print('',file=fout)
        with open(f'ROMAN_NIR_cadence{cadence:.0f}.IGNORE','w') as fout:
            print('',file=fout)
        os.chdir('../')
    return

def edit_lcs_phase():
    # function for making the lightcurve versions with phase
    cadence = 3
    for minphase in [-5,0,5,10,15,20,25]:
        if not os.path.exists(f"ROMAN_NIR_phase{str(int(minphase)).replace('-','m')}"):
            os.system(f"mkdir ROMAN_NIR_phase{str(int(minphase)).replace('-','m')}")
        lcfiles = glob.glob('/usr/local/SNDATA_ROOT/SCRATCH/ROMAN_NIR/ROMAN_NIR_SN*.DAT')

        for l in lcfiles:
            #sn = snana.SuperNova(l)

            sn = snana.SuperNova(l.replace('/usr/local/SNDATA_ROOT/SCRATCH/ROMAN_NIR',
                                           f"ROMAN_NIR_phase{str(int(minphase)).replace('-','m')}"))
            sn.FLT = sn.BAND
            nobs = len(sn.FLT)
            #nobs = 0
            with open(l) as fin, \
                 open(l.replace('/usr/local/SNDATA_ROOT/SCRATCH/ROMAN_NIR',
                                f"ROMAN_NIR_phase{str(int(minphase)).replace('-','m')}"),'w') as fout:
                lastmjdZ,lastmjdY,lastmjdJ,lastmjdH = None,None,None,None
                for line in fin:
                    if line.startswith('NOBS:'):
                        print(f'NOBS: {nobs:.0f}',file=fout)
                    elif not line.startswith('OBS:'):
                        print(line.replace('\n',''),file=fout)
                    else:
                        mjd = float(line.split()[1])
                        #import pdb; pdb.set_trace()
                        if (mjd-float(sn.PEAKMJD.split()[0]))/(1+sn.z) < minphase: continue
                        filt = line.split()[2]
                        if filt == 'Z': lastmjd = lastmjdZ
                        elif filt == 'Y': lastmjd = lastmjdY
                        elif filt == 'J': lastmjd = lastmjdJ
                        elif filt == 'H': lastmjd = lastmjdH

                        if filt == 'Z' and lastmjd is None:
                            lastmjdZ = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'Z' and mjd >= lastmjdZ+cadence:
                            lastmjdZ = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'Y' and lastmjd is None:
                            lastmjdY = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'Y' and mjd >= lastmjdY+cadence:
                            lastmjdY = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'J' and lastmjd is None:
                            lastmjdJ = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'J' and mjd >= lastmjdJ+cadence:
                            lastmjdJ = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'H' and lastmjd is None:
                            lastmjdH = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        elif filt == 'H' and mjd >= lastmjdH+cadence:
                            lastmjdH = float(line.split()[1])
                            print(line.replace('\n',''),file=fout)
                        
        os.chdir(f"ROMAN_NIR_phase{str(int(minphase)).replace('-','m')}")
        os.system(f"ls * > ROMAN_NIR_phase{str(int(minphase)).replace('-','m')}.LIST")
        with open(f"ROMAN_NIR_phase{str(int(minphase)).replace('-','m')}.README",'w') as fout:
            print('',file=fout)
        with open(f"ROMAN_NIR_phase{str(int(minphase)).replace('-','m')}.IGNORE",'w') as fout:
            print('',file=fout)
        os.chdir('../')
    return


def mkcuts(fr):

    iCut = (fr.AV < 1) & (fr.STRETCH > 0.8) & (fr.STRETCH < 1.3) & (fr.STRETCHERR < 0.1) & (fr.PKMJDERR < 2*(1+fr.zHD)) &\
           (fr.AVERR < 0.1) & (fr.AV > -1) & (fr.FITPROB > 0.01)
    for k in fr.__dict__.keys():
        fr.__dict__[k] = fr.__dict__[k][iCut]
    return fr

def main(sigint=0.05):

    plt.subplots_adjust(wspace=0.5)

    ax = plt.subplot(211)
    ax.set_xlabel('cadence (days)',fontsize=13)
    ax.set_ylabel('RMS (mag)',fontsize=13)

    rmslist = []
    cadencelist = [1,3,5,7,9,11,13]
    for cadence in cadencelist:
        fr = txtobj(f'fitres/ROMAN_NIR_cadence{cadence:.0f}.FITRES.TEXT',fitresheader=True)
        #print(len(fr.CID))
        fr = mkcuts(fr)
        #print(len(fr.CID))
        mures = fr.DLMAG - cosmo.mu(fr.zHD)
        #import pdb; pdb.set_trace()
        rmslist += [np.sqrt(np.std(mures)**2.+sigint**2.)]
        rmserr = np.std([np.std(np.random.choice(mures,size=len(mures))) for i in range(100)])
        #rmslist += [np.std(mures)]
        #break
        #import pdb; pdb.set_trace()
    ax.errorbar(cadencelist,rmslist,yerr=rmserr,fmt='o-',color='k')
    axright = ax.twinx()
    axright.yaxis.tick_right()
    left_ylims = np.array(ax.get_ylim())
    axright.set_ylim((left_ylims/rmslist[0]-1.0)*100)
    axright.set_ylabel('% increase',fontsize=13)
    
    ax2 = plt.subplot(212)
    rmslist = []
    minphaselist = [-5,0,5,10,15,20,25]
    for minphase in minphaselist:
        fr = txtobj(f"fitres/ROMAN_NIR_phase{str(int(minphase)).replace('-','m')}.FITRES.TEXT",fitresheader=True)
        fr = mkcuts(fr)
        mures = fr.DLMAG - cosmo.mu(fr.zHD)
        rmslist += [np.sqrt(np.std(mures)**2.+sigint**2.)]
        rmserr = np.std([np.std(np.random.choice(mures,size=len(mures))) for i in range(100)])
        #import pdb; pdb.set_trace()
    ax2.errorbar(minphaselist,rmslist,yerr=rmserr,fmt='o-',color='k')
    ax2.set_xlabel('min. phase (days)',fontsize=13)
    ax2.set_ylabel('RMS (mag)',fontsize=13)
    ax2right = ax2.twinx()
    ax2right.yaxis.tick_right()
    left_ylims = np.array(ax2.get_ylim())
    ax2right.set_ylim((left_ylims/rmslist[0]-1.0)*100)
    ax2right.set_ylabel('% increase',fontsize=13)

    #for axtmp in [ax,ax2]:
    #    axtmp.tick_params(top="on",bottom="on",left="on",right="off",direction="inout",length=8, width=1.5)
    #for axtmp in [axright,ax2right]:
    #    axtmp.tick_params(top="on",bottom="on",left="off",right="on",direction="inout",length=8, width=1.5)
    
    import pdb; pdb.set_trace()

def run_lcfit():
    os.system('snlc_fit.exe nml/snfit_ROMAN_CADENCE1.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_CADENCE3.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_CADENCE5.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_CADENCE7.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_CADENCE9.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_CADENCE11.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_CADENCE13.nml')

def run_lcfit_phase():
    os.system('snlc_fit.exe nml/snfit_ROMAN_PHASEm5.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_PHASE0.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_PHASE5.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_PHASE5.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_PHASE10.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_PHASE15.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_PHASE20.nml')
    os.system('snlc_fit.exe nml/snfit_ROMAN_PHASE25.nml')
    
if __name__ == "__main__":
    main()
    #edit_simlib()
    #edit_lcs()
    #edit_lcs_phase()
    #run_lcfit()
    #run_lcfit_phase()
