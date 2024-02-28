Implementations and parameter selection of Plover, a Masking-friendly hash-and-sign lattice
signature
============================================

This repository includes C and Python implementations of Plover-RLWE and Plover-NTRU,
along with Python scripts in the folder `params` to perform their parameter selection.

Plover was introduced in

> Muhammed Esgin, Thomas Espitau, Guilhem Niot, Thomas Prest, Amin
Sakzad, and Ron Steinfeld. Plover: Masking-friendly hash-and-sign lattice
signatures. In EUROCRYPT, 2024.

Note that a large part of the codebase is reused from [the
reference implementations of the NIST Additional Signature candidate Raccoon](https://github.com/masksign/raccoon).