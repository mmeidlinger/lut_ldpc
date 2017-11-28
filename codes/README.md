# Usage
Inventory of LDPC codes stored in  `.alist` format. Any parity check matrix specified in some file 
`filename.alist` migh also have its corresponding generator saved as `filename.gen.bin`.

# List of codes
Please update this file if you add new codes to the Repo.

| Filename                     | Description  | Source
|------------------------------|:-------------|:----
| [rate0.5_irreg_dvbs2_N64800.alist](rate0.5_irreg_dvbs2_N64800.alist)     |  Irregular LDPC Code used in the DVB-S2 Standarf | Extracted from MATLAB using [matlab_matrix_to_alist](../../sourcecode/decoder/matlab_matrix_to_alist.cpp)
| [rate0.84_reg_v6c32_N2048.alist](rate0.84_reg_v6c32_N2048.alist) |Regular (6,32) Code used in the Gigabit Ethernet Standard Note: The parity check of this code is rank deficient, so the rate is actually larger than 1-6/32| Gigabit Ethernet Standard.
| [rate0.84_irreg_v2-6c32_N2048.alist](rate0.84_irreg_v3-6c32_N2048.alist) | Irregular code, Gigabit ethernet standard code without linearly dependent rows | [alist2fullrank](../prog/alist2fullrank.cpp)  on [rate0.84_reg_v6c32_N2048.alist](rate0.84_reg_v6c32_N2048.alist) with seed 10
| [rate0.50_dv02-17_dc08-09_lut_q4_N10000.alist](rate0.50_dv02-17_dc08-09_lut_q4_N10000.alist) | Irregular Code optimized for LUTs with 4 bits resolution | [peg](../scripts/peg.sh) on ensemble [rate0.50_dv02-17_dc08-09_lut_q4.ens](../ensembles/rate0.50_dv02-17_dc08-09_lut_q4.ens)
| [rate0.50_dv02-17_dc08-09_lut_q4_N500.alist](rate0.50_dv02-17_dc08-09_lut_q4_N500.alist) | Short irregular Code for testing purposes | [peg](../scripts/peg.sh) on ensemble [rate0.50_dv02-17_dc08-09_lut_q4.ens](../ensembles/rate0.50_dv02-17_dc08-09_lut_q4.ens)
| [rate0.50_dv03_dc06_N10000.alist](rate0.50_dv03_dc06_N10000.alist) | (3,6) Regular Code | [peg](../scripts/peg.sh) on ensemble [rate0.50_dv03_dc06.ens](../ensembles/rate0.50_dv03_dc06.ens)rate0.50_dv02-17_dc08-09_lut_q4_PEG_N10000.alist
