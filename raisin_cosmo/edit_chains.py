#!/usr/bin/env python
# D. Jones - 7/20/21
# edit the planck cosmosis chains to add a varying SN magnitude param that
# has no impact on the likelihoods.

import numpy as np
planck_file_in = "cosmosis/planck_samples/chain_p-TTTEEE-lowE_wcdm.txt"
planck_file_out = "cosmosis/planck_samples/chain_p-TTTEEE-lowE_SNmag_wcdm.txt"

colnames = np.array(["#cosmological_parameters--omega_m","cosmological_parameters--h0","cosmological_parameters--omega_b","cosmological_parameters--n_s","cosmological_parameters--a_s","cosmological_parameters--omnuh2","cosmological_parameters--w","cosmological_parameters--tau","supernova_params--m","planck--a_planck","planck--a_cib_217","planck--xi_sz_cib","planck--a_sz","planck--ps_a_100_100","planck--ps_a_143_143","planck--ps_a_143_217","planck--ps_a_217_217","planck--ksz_norm","planck--gal545_a_100","planck--gal545_a_143","planck--gal545_a_143_217","planck--gal545_a_217","planck--calib_100t","planck--calib_217t","planck--galf_te_a_100","planck--galf_te_a_100_143","planck--galf_te_a_100_217","planck--galf_te_a_143","planck--galf_te_a_143_217","planck--galf_te_a_217","COSMOLOGICAL_PARAMETERS--SIGMA_8","COSMOLOGICAL_PARAMETERS--SIGMA_12","DATA_VECTOR--2PT_CHI2","prior","like","post","weight"])
iSN = np.where(colnames == "supernova_params--m")[0]

def main():

    with open(planck_file_in) as fin, open(planck_file_out,"w") as fout:
        for line in fin:
            line = line.replace("\n","")
            if line.startswith("#cosmological_parameters--omega_m"):
                print("\t".join(colnames.tolist()),file=fout)
            elif line.startswith("#"):
                print(line,file=fout)
            else:
                lineparts = line.split()
                newline = []
                for i,c in enumerate(colnames):
                    if i < iSN: newline += [lineparts[i]]
                    elif i == iSN: newline += ['%.13f'%np.random.uniform(-19.5,-18.5,1)]
                    elif i > iSN: newline += [lineparts[i-1]]
                print("\t".join(newline),file=fout)

if __name__ == "__main__":
    main()
