"""
test_ntt.py
Copyright (c) 2023 Plover Signature Team. See LICENSE.

=== Code for re-creating the NTT magic constants. Unit tests for NTT.
"""

from random import randrange
from polyr import *

if (__name__ == "__main__"):
    def _rand_poly(n=PLOVERSIGN_N, q=PLOVERSIGN_Q):
        """(TESTING) Random polynomial."""
        return [ randrange(q) for _ in range(n) ]

    #   test convolutions
    for _ in range(5):
        f = _rand_poly()
        g = _rand_poly()

        f_ntt = ntt(f)
        g_ntt = ntt(g)

        # we want to compute all the inverse of g_i
        a = 1

        for gi in g_ntt:
            a = (a*gi)%PLOVERSIGN_Q

        ainv = pow(a, PLOVERSIGN_PHIQ-1, PLOVERSIGN_Q)
        print(ainv)

        aainv = (a*ainv)%PLOVERSIGN_Q # can be unmasked as either 1, or we reject anyway
        if aainv != 1: # one of the g_i is not invertible
            print("one of the gi is not inverible")
            continue
    
        h_ntt = [None]*PLOVERSIGN_N

        proda = [1]
        for gi in g_ntt:
            proda.append((proda[-1]*gi)%PLOVERSIGN_Q)

        prodarev = 1
        for i in reversed(range(PLOVERSIGN_N)):
            # compute the inverse of g_i
            gi_inv = (ainv*prodarev*proda[i])%PLOVERSIGN_Q
            h_ntt[i] = (f_ntt[i]*gi_inv)%PLOVERSIGN_Q

            prodarev = (prodarev*g_ntt[i])%PLOVERSIGN_Q

        print(h_ntt)

        # for gi in g_ntt:
        #     gi = PLOVERSIGN_Q1
        #     print((pow(gi, PLOVERSIGN_PHIQ-1, PLOVERSIGN_Q))%PLOVERSIGN_Q)
