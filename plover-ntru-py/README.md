#   Plover/ntrusign-py

Copyright (c) 2023 Plover Signature Team. See LICENSE.
*(Code adapted from Python reference code of Raccoon, originally written by Thomas Prest and Markku-Juhani O. Saarinen.)*

Python implementation of Masked Plover, aimed at readability.

Test: `python3 test_genkat.py`: Generate .rsp KAT files for all variants.

```
Makefile            Only target: make clean
mask_random.py      Dummy Mask Random Generator (LFSR127 MRG)
nist_kat_drbg.py    NIST KAT Generator DRBG
polyr.py            Polynomial ring arithmetic + NTT code
plover_api.py          Serializaton/deserialization, NIST style functions
plover_core.py         Plover signature scheme -- Core Algorithm.
README.md           This file
requirements.txt    Python requirements: pycryptodome
test_genkat.py      Generate NIST KAT files (.rsp files)
test_ntt.py         Basic tests for NTT, also generate tweak tables
```

