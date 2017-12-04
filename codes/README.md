# Codes

Inventory of LDPC codes stored in  `.alist` format. Any parity check matrix specified in some file 
`filename.alist` migh also have its corresponding generator saved as `filename.gen.bin`.

# List of codes
Please update this file if you add new codes to the Repo.

| Filename                     | Description
|--------------------------|:------------
| [rate0.5_irreg_dvbs2_N64800.alist](rate0.84_reg_v6c32_N2048.alist)     |  Irregular LDPC Code used in the DVB-S2 Standard
| [rate0.84_reg_v6c32_N2048.alist](rate0.84_reg_v6c32_N2048.alist) | Regular (6,32) Code used in the Gigabit Ethernet Standard Note: The parity check of this code is rank deficient, so the rate is actually larger than 1-6/32
| [rate0.50_dv02-15_dc08-09_N10000.alist](rate0.50_dv02-15_dc08-09_N10000.alist) | Irregular Code optimized for BP decoding
| [rate0.50_dv02-08_dc07-08_lut_q4_N64800.alist](rate0.50_dv02-08_dc07-08_lut_q4_N64800.alist) | Irregular Code optimized for LUTs with 4 bits resolution
| [rate0.50_dv02-17_dc08-09_lut_q4_N10000.alist](rate0.50_dv02-17_dc08-09_lut_q4_N10000.alist) | Irregular Code optimized for LUTs with 4 bits resolution
| [rate0.50_dv02-17_dc08-09_lut_q4_N1000.alist](rate0.50_dv02-17_dc08-09_lut_q4_N1000.alist) | Irregular Code optimized for LUTs with 4 bits resolution
| [rate0.50_dv02-17_dc08-09_lut_q4_N500.alist](rate0.50_dv02-17_dc08-09_lut_q4_N500.alist) | Short irregular Code for testing purposes
| [rate0.50_dv03_dc06_N1000.alist](rate0.50_dv03_dc06_N1000.alist) | (3,6) Regular Code
| [rate0.50_dv03_dc06_N10000.alist](rate0.50_dv03_dc06_N10000.alist) | (3,6) Regular Code
