# Ensembles
Inventory of LDPC ensembles stored in  `.ens` format.

# File Format

```
dv_act dc_act
variable_node_degrees
variable_node_edge_pmf
check_node_degrees
check_node_edge_pmf
```
The first line specifies the number of active degrees (non-zero entries of the degree distributions)

# List of codes
Please update this file if you add new ensembles to the Repo.

| Filename                     | Description
|:-----------------------------|:--------------------------------------------------------
| [rate0.50_dv02-04_dc05-06.ens](rate0.50_dv02-04_dc05-06.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree  4, cf. [[1]](#literature)
| [rate0.50_dv02-05_dc06-07.ens](rate0.50_dv02-05_dc06-07.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree  5, cf. [[1]](#literature)
| [rate0.50_dv02-08_dc06-07.ens](rate0.50_dv02-08_dc06-07.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree  8, cf. [[1]](#literature)
| [rate0.50_dv02-11_dc07-08.ens](rate0.50_dv02-11_dc07-08.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree 11, cf. [[1]](#literature)
| [rate0.50_dv02-15_dc08-09.ens](rate0.50_dv02-15_dc08-09.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree 15, cf. [[1]](#literature)
| [rate0.50_dv02-50_dc09-11.ens](rate0.50_dv02-50_dc09-11.ens) | Irregular, rate 0.5 ensemble maximized for BIAWGN channel and maximum VN degree 50, cf. [[1]](#literature)
| [rate0.50_dv02-17_dc08-09_lut_q4.ens](rate0.50_dv02-17_dc08-09_lut_q4.ens) |  Irregular, rate 0.5 ensemble maximized for BIAWGN channel folowed by a 4 bit LUT decoder, cf. [[2]](#literature)
| [rate0.50_dv02-08_dc07-08_lut_q4.ens](rate0.50_dv02-08_dc07-08_lut_q4.ens) |  Irregular, rate 0.5 ensemble maximized for BIAWGN channel folowed by a 4 bit LUT decoder, cf. [[2]](#literature)
| [rate0.50_dv03_dc06.ens](rate0.50_dv03_dc06.ens) | Regular, rate 0.5 (3,6) ensemble
|[rate0.84_dv02-06_dc32.ens](rate0.84_dv02-06_dc32.ens)| Irregular, rate 0.84 ensemble


# Literature
[[1] T. Richardson, M. Shokrollahi, and R. Urbanke, “Design of capacity-approaching irregular low-density parity-check codes,” IEEE Trans. Information Theory, vol. 47, no. 2, pp. 619–637, Feb. 2001.](
    http://ieeexplore.ieee.org/document/910578/?arnumber=910578)
    
[[2] M. Meidlinger and G. Matz, “On irregular LDPC codes with quantized message passing decoding,” in Proc. IEEE SPAWC 2017, Sapporo, Japan, Jul. 2017.
    ](http://ieeexplore.ieee.org/document/8227780/)

