This folder contains scripts to define concrete parameters of Plover-RLWE
and Plover-NTRU.

It leverages the [lattice estimator](https://github.com/malb/lattice-estimator)
to estimate hardness of MLWE. Overstretched NTRU hardness is estimated with parameters obtained from https://github.com/WvanWoerden/NTRUFatigue.

Scripts `param-ring.py` and `param-ntru.py` compute parameters for unmasked
Plover, and the corresponding RLWE, RSIS and NTRU hardness.

Scripts `concrete-param-ntru.py` and `concrete-param-ring.py` compute parameters
for masked Plover, and concrete values for the modulus `q` enabling NTT and
computations over 32-bit integers, and concrete width for the uniform
distributions for masking orders `d` among 1, 2, 4, 8, 16, 32. They also 
generate header files to integrate parameters in our reference implementations.