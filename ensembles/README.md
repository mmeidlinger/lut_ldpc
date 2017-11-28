# Usage
Inventory of LDPC ensembles stored in  `.ens` format. 

# File Format
```
ds_act dvs_act dvc_act dc_act
symbol_node_degrees
symbol_node_edge_pmf
variable_to_symbol_node_degrees
variable_to_symbol_node_edge_pmf
variable_to_check_node_degrees
variable_to_check_node_edge_pmf
check_node_degrees
check_node_edge_pmf
```
The first line specifies the number of active degrees (non-zero entries of the degree distributions)
For conventional LDPC codes, 
`ds_act = dvs_act = symbol_node_degrees = variable_to_symbol_node_degrees  = 1` 
and 
`symbol_node_edge_pmf = variable_to_symbol_node_edge_pmf = 1.0`

# List of codes
Please update this file if you add new ensembles to the Repo.

| Filename                     | Description                                              | Source
|:-----------------------------|:---------------------------------------------------------|:---------
| [rate0.50_dv02-04_dc05-06.ens](rate0.50_dv02-04_dc05-06.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree  4 | 1 
| [rate0.50_dv02-05_dc06-07.ens](rate0.50_dv02-05_dc06-07.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree  5 | 1
| [rate0.50_dv02-09_dc06-08.ens](rate0.50_dv02-09_dc06-08.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree  9 | 1
| [rate0.50_dv02-11_dc07-08.ens](rate0.50_dv02-11_dc07-08.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree 11 | 1
| [rate0.50_dv02-15_dc08-09.ens](rate0.50_dv02-15_dc08-09.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree 15 | 1
| [rate0.50_dv02-50_dc09-11.ens](rate0.50_dv02-50_dc09-11.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree 50 | 1
| [rate0.50_dv02-17_dc08-09_lut_q4.ens](rate0.50_dv02-17_dc08-09_lut_q4.ens) |  Irregular, rate 0.5 ensemble maximized for BIAWGN channel folowed by a 4 bit LUT decoder | [degree_opt](../prog/degree_opt.cpp)
| [rate0.5_irreg_dvbs2.ens](rate0.5_irreg_dvbs2.ens)| Irregular, rate 0.5 ensemble of dvbs2 standard, removed degree 1 mass and renormalized  | [gen_ensemble](../prog/gen_ensemble.cpp) on [rate0.5_irreg_dvbs2_N64800.alist](../codes/rate0.5_irreg_dvbs2_N64800.alist)
| [rate0.50_dv03_dc06.ens](rate0.50_dv03_dc06.ens) | Regular, rate 0.5 (3,6) ensemble | Trivial
| [rate0.50_dv49_dc98.ens](rate0.50_dv49_dc98.ens) | Regular, rate 0.5 (49,98) ensemble. This was to test wether LUT decoders get worse for high degrees | Trivial
|[rate0.84_dv02-06_dc32.ens](rate0.84_dv02-06_dc32.ens)| Irregular, rate 0.84 ensemble| [gen_ensemble](../prog/gen_ensemble.cpp) on [rate0.84_irreg_v2-6c32_N2048.alist](rate0.84_irreg_v3-6c32_N2048.alist)


## Literature
1. [Design of capacity-approaching irregular low-density parity-check codes, by Richardson et al., TIT2001](
    http://ieeexplore.ieee.org/document/910578/?arnumber=910578)

