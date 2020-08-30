#!/usr/bin/env python
import numpy as np
import pylab as plt
import f90nml
import argparse
import os

_nir_nml = ['$RAISIN_ROOT/cosmo/fit/CSP_RAISIN.nml',
            '$RAISIN_ROOT/cosmo/fit/PS1_RAISIN.nml',
            '$RAISIN_ROOT/cosmo/fit/DES_RAISIN.nml']
_outdirs = ['$RAISIN_ROOT/cosmo/output/fit_nir_sys/CSP_RAISIN',
            '$RAISIN_ROOT/cosmo/output/fit_nir_sys/PS1_RAISIN',
            '$RAISIN_ROOT/cosmo/output/fit_nir_sys/DES_RAISIN']

_fitopt_dict = {'MWEBV':('MWEBV_SCALE 0.95','MWEBV_SCALE 0.95','MWEBV_SCALE 0.95'),
                'HST_CAL':('MAGOBS_SHIFT_ZP_PARAMS 0 0.00714 0',
                           'MAGOBS_SHIFT_ZP_PARAMS 0 0.00714 0',
                           'MAGOBS_SHIFT_ZP_PARAMS 0 0.00714 0'),
                'VPEC':('VPEC_FILE \'$RAISIN_ROOT/cosmo/vpec_sys_raisin.list\'',
                        'VPEC_FILE \'$RAISIN_ROOT/cosmo/vpec_sys_raisin.list\'',
                        'VPEC_FILE \'$RAISIN_ROOT/cosmo/vpec_sys_raisin.list\''),
                'MASS_DIVIDE':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'CSP_Y_SURVCAL':('MAGOBS_SHIFT_ZP \'Y 0.01 y 0.01\'','->FITOPT000','->FITOPT000'),
                'CSP_J_SURVCAL':('MAGOBS_SHIFT_ZP \'J 0.01 j 0.01\'','->FITOPT000','->FITOPT000'),
                'CSP_H_SURVCAL':('MAGOBS_SHIFT_ZP \'H 0.01 h 0.01\'','->FITOPT000','->FITOPT000'),
                'BIASCOR_SHAPE':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'BIASCOR_AV':('->FITOPT000','->FITOPT000','->FITOPT000'),
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

        return parser
        
    def mk_nml(self):

        for i,nml in enumerate(_nir_nml):
            nml = os.path.expandvars(nml)
            with open(nml.replace('.nml','_sys.nml'),'w') as fout:
                print(f'OUTDIR: {_outdirs[i]}',file=fout)
                with open(nml) as fin:
                    for line in fin:
                        if not line.startswith('OUTDIR'): print(line.replace('\n',''),file=fout)
                print('',file=fout)
                for k in _fitopt_dict.keys():
                    # no kcor variants for now
                    if 'KCOR' in k: continue
                    
                    print(f'FITOPT: [{k}] {_fitopt_dict[k][i]}',file=fout)

    def mk_cosmo_inputs(self):
        raise NotImplementedError('')
                    
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
