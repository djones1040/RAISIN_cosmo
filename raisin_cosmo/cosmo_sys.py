#!/usr/bin/env python
import numpy as np
import pylab as plt
import f90nml
import argparse
import os
import glob
import snana
from raisin_cosmo import get_vpec
from scipy.optimize import minimize
from txtobj import txtobj
import copy
import cosmo
import scipy.stats
from scipy.stats import binned_statistic
from raisin_cosmo import ovdatamc
#export PYTHONPATH=$SNANA_DIR/util:$PYTHONPATH

_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_CIDS.LIST',unpack=True,dtype=str),
                            #np.loadtxt('output/goodcids/CfA_CIDS.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/PS1_CIDS.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/DES_CIDS.LIST',unpack=True,dtype=str)))

_cosmomc_batch = """#!/bin/bash
#SBATCH --job-name=RAISIN_stat
#SBATCH --time=34:00:00
###SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=1-2
#SBATCH --cpus-per-task=1
#SBATCH --partition=broadwl-lc
#SBATCH --output=/scratch/midway2/rkessler/djones/cosmomc/chains/RAISIN_stat.log
#SBATCH --account=pi-rkessler
#SBATCH --mem=20GB

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module unload openmpi
module load intelmpi/5.1+intel-16.0
module load cfitsio/3
module load mkl
module load gcc/6.1
source /project2/rkessler/SURVEYS/PS1MD/USERS/djones/RAISIN/code/plc_3.0/plc-3.01/bin/clik_profile.sh

PARAMS=`expr ${{SLURM_ARRAY_TASK_ID}} - 1`

INI_FILES=({} {})
DONE_FILES=(done_0.txt done_1.txt)

cd /scratch/midway2/rkessler/djones/cosmomc/chains/
#mpirun /project2/rkessler/PRODUCTS/CosmoMC/v03/CosmoMC-master/cosmomc ${{INI_FILES[$PARAMS]}}
mpirun $RAISIN_ROOT/CosmoMC/cosmomc ${{INI_FILES[$PARAMS]}}

if [ $? -eq 0 ]; then
    echo "SUCCESS" > ${{DONE_FILES[$PARAMS]}}
else
    echo "FAILURE" > ${{DONE_FILES[$PARAMS]}}
fi

"""

_cosmomc_ini = """DEFAULT(batch2/JLA.ini)
DEFAULT(batch2/GAUSS.ini)
DEFAULT(batch2/common.ini)
accuracy_level = 1.35
high_accuracy_default = F
CMB_lensing = F
lmin_store_all_cmb = 0
stop_on_error = F

MPI_Converge_Stop = 0.002
MPI_Limit_Converge = 0.002
MPI_Limit_Converge_Err = 0.13
MPI_Max_R_ProposeUpdate = 0.03
indep_sample = 5

propose_matrix= planck_covmats/base_TT_lowTEB_plik.covmat
root_dir = chains/

action = 0
num_threads = 1

start_at_bestfit = F
feedback = 1
use_fast_slow = F
checkpoint = T

#VIN LETS TRY SIMPLE METROPOLIS WITH THAT LIKELIHOOD
#sampling_method = 1
#sampling_method=7 is a new fast-slow scheme good for Planck
sampling_method = 7
dragging_steps  = 3
propose_scale = 2

#these are just small speedups for testing
get_sigma8=F

param[wa]=0
param[w]=-0.995 -1.5 -0.5 0.001 0.001
param[omegak]=0
param[omegam]=0.3 0.2 0.4 0.002 0.002
compute_tensors=F
param[r]=0.0
param[calPlanck]=1
# VERY IMPORTANT TO FIX PARAMETERS WE DONT CONSTRAIN WITH OUR APPROXIMATION
param[tau]=0.078
param[logA]=3.090388
param[alpha_JLA]=0.14
param[beta_JLA]=3.1
file_root={}
jla_dataset={}"""

_datasettext = """#Settings for the joint SNLS/SDSS data analysis
name = JLA
data_file = {}
pecz = 0
intrinsicdisp = 0
intrinsicdisp0 = 0
intrinsicdisp1 = 0
intrinsicdisp2 = 0
intrinsicdisp3 = 0
twoscriptmfit = F
#scriptmcut = 10.0
has_mag_covmat = {}
mag_covmat_file = {}
has_stretch_covmat = F
has_colour_covmat = F
has_mag_stretch_covmat = F
has_mag_colour_covmat = F
has_stretch_colour_covmat = F"""

_nir_nml = ['$RAISIN_ROOT/cosmo/fit/CSP_RAISIN.nml',
            '$RAISIN_ROOT/cosmo/fit/PS1_RAISIN.nml',
            '$RAISIN_ROOT/cosmo/fit/DES_RAISIN.nml']
_outdirs = ['$RAISIN_ROOT/cosmo/output/fit_nir_sys/CSP_RAISIN',
            '$RAISIN_ROOT/cosmo/output/fit_nir_sys/PS1_RAISIN',
            '$RAISIN_ROOT/cosmo/output/fit_nir_sys/DES_RAISIN']
_data_dirs = ['$RAISIN_ROOT/cosmo/data/Photometry/CSPDR3_RAISIN',
              '$RAISIN_ROOT/cosmo/data/Photometry/PS1_RAISIN',
              '$RAISIN_ROOT/cosmo/data/Photometry/DES_RAISIN']

_sysgroupdict = {'photcal':('CSP_Y_SURVCAL','CSP_J_SURVCAL','CSP_H_SURVCAL','HST_CAL'),
                 'hstcal':('HST_CAL',),
                 'lowzcal':('CSP_Y_SURVCAL','CSP_J_SURVCAL','CSP_H_SURVCAL'),
                 'massdivide':('MASS_DIVIDE',),
                 'biascor':('BIASCOR_SHAPE_LOWZ','BIASCOR_AV_LOWZ','BIASCOR_SHAPE_HIGHZ','BIASCOR_AV_HIGHZ'),
                 'pecvel':('VPEC',),
                 'mwebv':('MWEBV',),
                 'kcor':('KCOR1','KCOR2','KCOR3','KCOR4','KCOR5','KCOR6','KCOR7','KCOR8','KCOR9','KCOR10')}


_fitopt_dict = {'MWEBV':('MWEBV_SCALE 0.95','MWEBV_SCALE 0.95','MWEBV_SCALE 0.95'),
                'HST_CAL':('MAGOBS_SHIFT_ZP_PARAMS 0 0.00714 0',
                           'MAGOBS_SHIFT_ZP_PARAMS 0 0.00714 0',
                           'MAGOBS_SHIFT_ZP_PARAMS 0 0.00714 0'),
                'VPEC':('HEADER_OVERRIDE_FILE \'$RAISIN_ROOT/cosmo/vpec_sys_raisin.list\'',
                        'HEADER_OVERRIDE_FILE \'$RAISIN_ROOT/cosmo/vpec_sys_raisin.list\'',
                        'HEADER_OVERRIDE_FILE \'$RAISIN_ROOT/cosmo/vpec_sys_raisin.list\''),
                'MASS_DIVIDE':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'MASS_STEP':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'CSP_Y_SURVCAL':('MAGOBS_SHIFT_ZP \'Y 0.01 y 0.01\'','->FITOPT000','->FITOPT000'),
                'CSP_J_SURVCAL':('MAGOBS_SHIFT_ZP \'J 0.01 j 0.01\'','->FITOPT000','->FITOPT000'),
                'CSP_H_SURVCAL':('MAGOBS_SHIFT_ZP \'H 0.01\'','->FITOPT000','->FITOPT000'),
                'BIASCOR_SHAPE_LOWZ':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'BIASCOR_AV_LOWZ':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'BIASCOR_SHAPE_HIGHZ':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'BIASCOR_AV_HIGHZ':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'KCOR1':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_KCOR1.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_KCOR1.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_KCOR1.fits\''),
                'KCOR2':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_KCOR2.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_KCOR2.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_KCOR2.fits\''),
                'KCOR3':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_KCOR3.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_KCOR3.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_KCOR3.fits\''),
                'KCOR4':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_KCOR4.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_KCOR4.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_KCOR4.fits\''),
                'KCOR5':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_KCOR5.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_KCOR5.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_KCOR5.fits\''),
                'KCOR6':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_KCOR6.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_KCOR6.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_KCOR6.fits\''),
                'KCOR7':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_KCOR7.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_KCOR7.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_KCOR7.fits\''),
                'KCOR8':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_KCOR8.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_KCOR8.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_KCOR8.fits\''),
                'KCOR9':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_KCOR9.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_KCOR9.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_KCOR9.fits\''),
                'KCOR10':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_KCOR10.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_KCOR10.fits\'',
                         'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_KCOR10.fits\''),
}

def writecov(covmat,covmatfile):

    fout = open(covmatfile,'w')
    print('%i'%np.shape(covmat)[0],file=fout)
    for i in range(np.shape(covmat)[0]):
        for j in range(np.shape(covmat)[1]):
            if i != j:
                print('%8.5e'%covmat[j,i],file=fout)
            else:
                print('%8.5e'%0,file=fout)
    fout.close()

class biascor:
    def __init__(self):

        self.nirsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN_NIR_SIM',
                                 '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN_NIR_SIM',
                                 '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN_NIR_SIM']
        self.nirsimfitreslist = [os.path.expandvars(filepath) for filepath in self.nirsimfitreslist]
        
        self.opticalnirsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_all/CSP_RAISIN_OPTNIR_SIM',
                                        '$RAISIN_ROOT/cosmo/output/fit_all/PS1_RAISIN_OPTNIR_SIM',
                                        '$RAISIN_ROOT/cosmo/output/fit_all/DES_RAISIN_OPTNIR_SIM']
        self.opticalnirsimfitreslist = [os.path.expandvars(filepath) for filepath in self.opticalnirsimfitreslist]


        self.nirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN.FITRES.TEXT',
                                  '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN.FITRES.TEXT',
                                  '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN.FITRES.TEXT']
        self.nirdatafitreslist = [os.path.expandvars(filepath) for filepath in self.nirdatafitreslist]

        self.opticalnirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',
                                         '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT',
                                         '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT']
        self.opticalnirdatafitreslist = [os.path.expandvars(filepath) for filepath in self.opticalnirdatafitreslist]

    def apply_all_cuts(self,fr,fropt,restrict_to_good_list=False):

        # AV
        iGoodAV = np.zeros(len(fr.CID),dtype=bool)
        for j,i in enumerate(fr.CID):
            if i in fropt.CID and fropt.AV[fropt.CID == i] < 0.3*fropt.RV[fropt.CID == i]:
                iGoodAV[j] = True

        # reasonable stretch
        iGoodSt = np.zeros(len(fr.CID),dtype=bool)
        for j,i in enumerate(fr.CID):
            if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > 0.8 and fropt.STRETCH[fropt.CID == i] < 1.175:
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

    def apply_biascor(self,fitopt,sys=None):
        if fitopt == 0:
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            fig2 = plt.figure()
            
        # can interpolate for each SN individually because samples are so small
        for nirdatadir,opticalnirdatafitres,nirsimfitres,opticalnirsimfitres,name,idx in zip(
                _outdirs,self.opticalnirdatafitreslist,
                self.nirsimfitreslist,self.opticalnirsimfitreslist,['CSP','PS1','DES'],range(4)):

            frdata = txtobj(glob.glob(os.path.expandvars("%s/*/FITOPT%03i.FITRES"%(nirdatadir,fitopt)))[0],fitresheader=True)
            fropt = txtobj(opticalnirdatafitres,fitresheader=True)
            with open(os.path.expandvars(f"{nirdatadir}/FITOPT.README")) as fin:
                fitoptstr = fin.readlines()[1:]
            sys = fitoptstr[fitopt]
            if 'BIASCOR_AV_LOWZ' in sys and 'CSP_RAISIN' in nirdatadir: syskey = '_AVSYS'; print(syskey,fitopt)
            elif 'BIASCOR_SHAPE_LOWZ' in sys and 'CSP_RAISIN' in nirdatadir: syskey = '_STSYS'; print(syskey,fitopt)
            elif 'BIASCOR_AV_HIGHZ' in sys and 'CSP_RAISIN' not in nirdatadir: syskey = '_AVSYS'; print(syskey,fitopt)
            elif 'BIASCOR_SHAPE_HIGHZ' in sys and 'CSP_RAISIN' not in nirdatadir: syskey = '_STSYS'; print(syskey,fitopt)
            else: syskey = ''

            frsim = txtobj(glob.glob(os.path.expandvars(f'{nirsimfitres}/*{syskey}/FITOPT000.FITRES'))[0],
                           fitresheader=True)
            froptsim = txtobj(glob.glob(os.path.expandvars(f'{opticalnirsimfitres}/*{syskey}/FITOPT000.FITRES'))[0],
                              fitresheader=True)
            
            frdata = self.apply_all_cuts(frdata,fropt,restrict_to_good_list=True)
            fropt = self.apply_all_cuts(fropt,fropt,restrict_to_good_list=True)
            frsim = self.apply_all_cuts(frsim,froptsim,restrict_to_good_list=False)
            froptsim = self.apply_all_cuts(froptsim,froptsim,restrict_to_good_list=False)
            
            frdata.DLMAG_biascor = np.array([-99.]*len(frdata.CID))
            for j,i in enumerate(frdata.CID):
                if name not in ['PS1','DES']: iBias = np.where(np.abs(frsim.zCMB - frdata.zCMB[j]) < 0.015)[0]
                else: iBias = np.where(np.abs(frsim.zCMB - frdata.zCMB[j]) < 0.03)[0]
                if len(iBias) < 60: #95:
                    import pdb; pdb.set_trace()
                    raise RuntimeError('not enough biascor events!')
                bias_j = np.average(frsim.DLMAG[iBias]-frsim.SIM_DLMAG[iBias],weights=1/(frsim.DLMAGERR[iBias]))
                frdata.DLMAG_biascor[j] = bias_j
                frdata.DLMAG[j] -= bias_j

            if fitopt == 0:
                zbins = np.linspace(np.min(frsim.zCMB),np.max(frsim.zCMB),10)
                biasbinned = binned_statistic(frsim.zCMB,frsim.DLMAG-frsim.SIM_DLMAG,bins=zbins,statistic='median').statistic
                ax.plot((zbins[1:]+zbins[:-1])/2.,biasbinned,'o-',label=name)
                fig2.clf()
                hist = ovdatamc.ovhist()
                parser = hist.add_options(usage='')
                options,  args = parser.parse_args()
                hist.options = options
                hist.options.journal = True
                hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.8,1.175),('AV',-0.5,0.93)]
                hist.options.nbins = 10
                hist.options.clobber = True
                hist.options.outfile = 'figs/sim_%s.png'%name
                hist.options.histvar = ['MURES','AV','STRETCH']
                hist.main(data=fropt,sim=froptsim)
                
            frdata.writefitres(f"output/fitres_cosmo/{name}.FITRES",clobber=True)
            if idx == 0:
                frdata_combined = copy.deepcopy(frdata)
            else:
                for k in frdata_combined.__dict__.keys():
                    frdata_combined.__dict__[k] = np.append(frdata_combined.__dict__[k],frdata.__dict__[k])

        frdata_combined.writefitres('output/cosmo_fitres/RAISIN_combined_FITOPT%03i.FITRES'%fitopt,clobber=True)
        frdata_combined.DLMAG += frdata_combined.DLMAG_biascor
        frdata_combined.writefitres('output/cosmo_fitres/RAISIN_combined_FITOPT%03i_nobiascor.FITRES'%fitopt,clobber=True)
        if fitopt == 0:
            ax.legend()
            fig.savefig('figs/biascor_applied.png')

class cosmo_sys:
    def __init__(self):
        pass
                
    def add_arguments(self,parser=None, usage=None, config=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler="resolve")

        # The basics
        parser.add_argument('-v', '--verbose', action="count", dest="verbose",
                            default=0,help='verbosity level')
        parser.add_argument('-n','--make_nml', default=False,action="store_true",
                            help='NML file with FITOPT')
        parser.add_argument('-c','--make_covmat', default=False, action="store_true",
                            help='covmat & cosmoMC inputs')
        parser.add_argument('--clobber', default=False,action="store_true",
                            help='clobber flag')

        return parser

    def get_vpec(self):
        with open(os.path.expandvars('$RAISIN_ROOT/cosmo/vpec_sys_raisin.list'),'w') as foutsys,\
             open(os.path.expandvars('$RAISIN_ROOT/cosmo/vpec_baseline_raisin.list'),'w') as foutbase:
            print("VARNAMES: CID VPEC VPEC_ERR",file=foutsys)
            print("VARNAMES: CID VPEC VPEC_ERR",file=foutbase)
            for d in _data_dirs:
                listfile = glob.glob(os.path.expandvars(f"{d}/*LIST"))[0]
                files = np.loadtxt(listfile,unpack=True,dtype='str')
                for f in files:
                    sn = snana.SuperNova(os.path.expandvars(f"{d}/{f}"))
                    if 'DECL' not in sn.__dict__.keys():
                        sn.DECL = sn.DEC
                    vpec,vpec_sys = get_vpec.main(float(sn.RA.split()[0]),float(sn.DECL.split()[0]),float(sn.REDSHIFT_FINAL.split('+-')[0]))
                    print(f"SN: {sn.SNID} {vpec+vpec_sys} 250",file=foutsys)
                    print(f"SN: {sn.SNID} {vpec} 250",file=foutbase)
                    
    def mk_nml(self):

        # peculiar velocity list
        if self.options.clobber or not os.path.exists(os.path.expandvars('$RAISIN_ROOT/cosmo/vpec_sys_raisin.list')) or not\
           os.path.exists(os.path.expandvars('$RAISIN_ROOT/cosmo/vpec_baseline_raisin.list')):
            self.get_vpec()
        
        for i,nml in enumerate(_nir_nml):
            nml = os.path.expandvars(nml)
            with open(nml.replace('.nml','_sys.nml'),'w') as fout:
                print(f'OUTDIR: {_outdirs[i]}',file=fout)
                with open(nml) as fin:
                    for line in fin:
                        if not line.startswith('OUTDIR') and not line.startswith('APPEND_FITRES'):
                            print(line.replace('\n',''),file=fout)
                print('',file=fout)
                for k in _fitopt_dict.keys():
                    # no kcor variants for now
                    if 'KCOR' in k: continue
                    
                    print(f'FITOPT: [{k}] {_fitopt_dict[k][i]}',file=fout)

    def bias_correct(self):
        bc = biascor()
        with open(os.path.expandvars(f"{_outdirs[0]}/FITOPT.README")) as fin:
            fitoptcount = 0
            for line in fin:
                if line.startswith('FITOPT'): fitoptcount += 1
        for i in range(fitoptcount):
            bc.apply_biascor(i)

    def nuisance_params(self,fitopt):

        # simultaneously fit for mass step, 3 sigint bins and 3 distance bins
        fr = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT%03i.FITRES'%fitopt,fitresheader=True)
        fr.mures = fr.DLMAG - cosmo.mu(fr.zCMB)
        fr.mures -= np.median(fr.mures)

        with open(os.path.expandvars(f"{_outdirs[0]}/FITOPT.README")) as fin:
            fitoptstr = fin.readlines()[1:]

        sys = fitoptstr[fitopt]
        if 'MASS_DIVIDE' in sys:
            msteploc = 10.15
        else:
            msteploc = 10
        p_lm = np.zeros(len(fr.CID))-99.
        for i in range(len(fr.CID)):
            if fr.HOST_LOGMASS_ERR[i] == 10: fr.HOST_LOGMASS_ERR[i] = 0.2
            if fr.HOST_LOGMASS[i] > msteploc:
                p_lm[i] = scipy.stats.norm.cdf(
                    msteploc,fr.HOST_LOGMASS[i],
                    fr.HOST_LOGMASS_ERR[i])*100.
            else:
                p_lm[i] = scipy.stats.norm.cdf(
                    msteploc,fr.HOST_LOGMASS[i],
                    fr.HOST_LOGMASS_ERR[i])*100.

        iCSP = fr.IDSURVEY == 5
        iPS1_midz = (fr.IDSURVEY == 15) & (fr.zCMB < 0.4306)
        iPS1_highz = (fr.IDSURVEY == 15) & (fr.zCMB >= 0.4306)
        iDES_midz = (fr.IDSURVEY == 10) & (fr.zCMB < 0.4306)
        iDES_highz = (fr.IDSURVEY == 10) & (fr.zCMB >= 0.4306)
        
        def neglnlikefunc(x,p_lm=None,mu_i=None,sigma_i=None):

            mu_lowz,mu_midz,mu_highz = x[0],x[1],x[2]
            sigint_csp,sigint_ps1,sigint_des = x[3],x[4],x[5]
            mass_step = x[6]
            
            # sigint split by sample, but distance split by redshift
            # each one with a low-mass and high-mass component
            loglike_csp = -np.sum(np.logaddexp(-(mu_i[iCSP]+mass_step-mu_lowz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                                               np.log((1-0.01*p_lm[iCSP])/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.))),
                                               -(mu_i[iCSP]-mu_lowz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                                               np.log((0.01*p_lm[iCSP])/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.)))))
            
            loglike_ps1_midz = -np.sum(np.logaddexp(-(mu_i[iPS1_midz]+mass_step-mu_midz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                                                    np.log((1-0.01*p_lm[iPS1_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.))),
                                               -(mu_i[iPS1_midz]-mu_midz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                                                    np.log((0.01*p_lm[iPS1_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.)))))
            
            loglike_ps1_highz = -np.sum(np.logaddexp(-(mu_i[iPS1_highz]+mass_step-mu_highz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                                     np.log((1-0.01*p_lm[iPS1_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.))),
                                               -(mu_i[iPS1_highz]-mu_highz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                                     np.log((0.01*p_lm[iPS1_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.)))))
            
            loglike_des_midz = -np.sum(np.logaddexp(-(mu_i[iDES_midz]+mass_step-mu_midz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                                                    np.log((1-0.01*p_lm[iDES_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.))),
                                                    -(mu_i[iDES_midz]-mu_midz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                                                    np.log((0.01*p_lm[iDES_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.)))))
            
            loglike_des_highz = -np.sum(np.logaddexp(-(mu_i[iDES_highz]+mass_step-mu_highz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                                     np.log((1-0.01*p_lm[iDES_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.))),
                                                     -(mu_i[iDES_highz]-mu_highz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                                     np.log((0.01*p_lm[iDES_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.)))))

            return loglike_csp + loglike_ps1_midz + loglike_ps1_highz + loglike_des_midz + loglike_des_highz
        
        md = minimize(neglnlikefunc,[0,0.01,0.02,0.09,0.1,0.11,0.1],args=(p_lm,fr.mures,fr.DLMAGERR))


        # apply sigint and host mass, along with lensing uncertainty
        # write the new fitres files
        fr.DLMAGERR[fr.IDSURVEY == 15] = np.sqrt(fr.DLMAGERR[fr.IDSURVEY == 15]**2. + md.x[4]**2.)
        fr.DLMAGERR[fr.IDSURVEY == 5] = np.sqrt(fr.DLMAGERR[fr.IDSURVEY == 5]**2. + md.x[3]**2.)
        fr.DLMAGERR[fr.IDSURVEY == 10] = np.sqrt(fr.DLMAGERR[fr.IDSURVEY == 10]**2. + md.x[5]**2.)

        if 'MASS_STEP' not in sys:
            fr.DLMAG[fr.HOST_LOGMASS > 10] += md.x[6]/2.
            fr.DLMAG[fr.HOST_LOGMASS <= 10] -= md.x[6]/2.
        else:
            fr.DLMAG[fr.HOST_LOGMASS > 10] += (md.x[6]+np.sqrt(md.hess_inv[6,6]))/2.
            fr.DLMAG[fr.HOST_LOGMASS <= 10] -= (md.x[6]+np.sqrt(md.hess_inv[6,6]))/2.

        fr.writefitres('output/cosmo_fitres/RAISIN_combined_FITOPT%03i_new.FITRES'%fitopt,clobber=True) 
        #if 'MASS_DIVIDE' in sys: import pdb; pdb.set_trace()
        
    def sys_covmat(self):
        syslist = ['stat','all','photcal','hstcal','lowzcal',
                   #'massstep','kcor',
                   'massdivide',
                   'biascor','pecvel','mwebv']
        for sys in syslist:
            count = 0
            fin = open(os.path.expandvars(f'{_outdirs[0]}/FITOPT.README'),'r')
            syslines = fin.readlines()[1:-1]
            fin.close()

            for sysline,i in zip(syslines,range(len(syslines))):
                sysname = sysline.split('[')[-1].split(']')[0]
                if sys == 'all' or sys == 'stat' or sysname in _sysgroupdict[sys]:
                    if count == 0:
                        baselinefitres = 'output/cosmo_fitres/RAISIN_combined_FITOPT%03i_new.FITRES'%0
                        frbase = txtobj(baselinefitres,fitresheader=True)
                        covshape = len(frbase.CID)
                        
                        basecov = np.zeros([covshape,covshape])
                        for j in range(covshape):
                            basecov[j,j] = frbase.DLMAGERR[j]**2.

                        outcov = basecov[:]
                        #if sys == 'massdivide': import pdb; pdb.set_trace()
                        if sys == 'stat': break
                    count += 1
                    
                    fr = txtobj('output/cosmo_fitres/RAISIN_combined_FITOPT%03i_new.FITRES'%i,fitresheader=True)

                    global_off = np.average(fr.DLMAG-frbase.DLMAG,weights=1/(fr.DLMAGERR**2.))
                    dm2=fr.DLMAG-global_off-frbase.DLMAG
                    #if "/SALT2/" in sysline: dm2 *= 0.3
                    dm2t=np.matrix(dm2)
                    dm2t=dm2t.T
                    dmm=dm2t*np.matrix(dm2)
                    outcov = np.add(outcov,dmm)
                    #if sys == 'biascor': import pdb; pdb.set_trace()
                    #if 'MASS_DIVIDE' in sysline and sys == 'all': import pdb; pdb.set_trace()
                #outcov = np.add(basecov,outcov)

            writecov(outcov,'output/cosmo_fitres/RAISIN_%s.covmat'%sys)

            with open('output/cosmo_fitres/RAISIN_%s_lcparams.txt'%sys,'w') as fout:
                print('# name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec biascor snana',file=fout)
                for i in range(covshape):
                    print('0.0 %.6f %.6f 0.000000 %.6f %.6f 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.0000e+00 0.0000e+00 0.0000e+00 0.0000 0.0000000 0.0000000 0.0000 0'%(
                        frbase.zCMB[i],frbase.zCMB[i],frbase.DLMAG[i],np.sqrt(outcov[i,i])),file=fout)
                fout.close()

    def cosmomc_inputs(self,dosubmit=False):
        syslist = ['stat','all','photcal','hstcal','lowzcal',
                   'massstep','massdivide',
                   'biascor','pecvel','mwebv','kcor']

        for sys in syslist:
            batchfile = f'cosmomc/RAISIN_{sys}.sbatch'
            inifile = f'cosmomc/RAISIN_{sys}.ini'
            datasetfile = f'cosmomc/RAISIN_{sys}.dataset'
            lcparfile = f'output/cosmo_fitres/RAISIN_{sys}_lcparams.txt'
            covfile = f'output/cosmo_fitres/RAISIN_{sys}.covmat'
            root = 'raisin_{sys}'
            
            # batch file
            with open(batchfile,'w') as fout:
                print(_cosmomc_batch.format(inifile,inifile),file=fout)
            
            # ini file
            with open(inifile,'w') as fout:
                print(_cosmomc_ini.format(root,datasetfile),file=fout)
            
            # dataset file
            has_covmat = 'T'
            with open(datasetfile,'w') as fout:
                print(_datasettext.format(lcparfile,has_covmat,covfile),file=fout)
            
            # lcparams and covmats should already be written

            if dosubmit:
                os.system(f'qsub {batchfile}')
        
    def mk_cosmo_inputs(self):

        # how many fit options?
        with open(os.path.expandvars(f"{_outdirs[0]}/FITOPT.README")) as fin:
            fitoptcount = 0
            for line in fin:
                if line.startswith('FITOPT'): fitoptcount += 1
                
        # make files for every FITOPT, bias-correct, compute sigint, apply mass step
        self.bias_correct()
        for i in range(fitoptcount):
            self.nuisance_params(i)
        
        # sys covmat from FITRES files
        # covmat for each individual systematic shift also
        self.sys_covmat()

        # cosmomc inputs
        self.cosmomc_inputs()
        
    def main(self):
        
        if self.options.make_nml:
            self.mk_nml()
        elif self.options.make_covmat:
            self.mk_cosmo_inputs()
        else:
            raise RuntimeWarning('doing nothing!')
        
if __name__ == "__main__":

    cs = cosmo_sys()
    parser = cs.add_arguments()
    args = parser.parse_args()
    cs.options = args
    cs.main()
