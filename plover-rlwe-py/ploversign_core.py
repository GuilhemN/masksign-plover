"""
plover_core.py
Copyright (c) 2023 Plover Signature Team. See LICENSE.

=== Masked Plover signature scheme: Core implementation.
"""

import os

from Crypto.Hash import SHAKE256
from nist_kat_drbg import NIST_KAT_DRBG
from mask_random import MaskRandom
from polyr import *
from math import ceil, sqrt, log, floor

BYTEORDER = "little"

class PloverSign:

    ### Public Interface

    #   initialize
    def __init__(self,  bitsec,
                        q, logdivide, nut, rep, ut, up, n, d,
                        masking_poly=MaskRandom().random_poly,
                        random_bytes=os.urandom, kappa=512):
        """Initialize a Plover instance."""

        self.name   =   f'Plover-{bitsec}-{d}'
        self.bitsec =   bitsec
        self.d      =   d
        self.q      =   q
        self.logdivide = logdivide
        self.q_bits =   self.q.bit_length()
        self.n      =   n
        self.nut    =   nut
        self.rep    =   rep
        self.ut     =   ut
        self.up     =   up

        self.sec    =   self.bitsec//8  # pre-image resistance, bytes
        self.crh    =   2*self.sec      # collision resistance, bytes
        self.as_sz  =   self.sec        # A seed size
        self.salt_sz  =   self.crh        # mu digest H(tr, m) size
        self.tr_sz  =   self.crh        # tr digest H(pk) size
        self.mk_sz  =   self.sec        # serialization "key" size

        self.masking_poly = masking_poly
        self.random_bytes = random_bytes

        #   calculate derived parmeters
        self._compute_metrics()

    def keygen(self):
        """Plover keypair generation."""

        #   --- 1.  seed <- {0,1}^kappa
        seed = self.random_bytes(self.as_sz)

        #   --- 2.  a := ExpandA(seed)
        a_ntt = ntt(self._expand_a(seed))

        #  --- 3.  [[(e, s)]] <- AddRepNoise([[(e, s)]], ut, rep)
        me = self._vec_add_rep_noise( self.ut, 0, self.rep )

        ms = self._vec_add_rep_noise( self.ut, 0, self.rep )
        ms_ntt = mat_ntt([ms])[0]

        #  --- 4.  [[b]] := A * [[s]] + [[e]]
        mb_ntt = mul_mat_mvec_ntt([[a_ntt]], [ms_ntt])[0]
        mb = mat_intt([mb_ntt])[0]
        for j in range(self.d):
            mb[j] = poly_add(mb[j], me[j])

        #   --- 7.  b := Unmask([[b]])
        b = self._decode(mb)

        #   --- 8.  b := round( eta - b )_q->q_t
        b[0] = (2**self.logdivide - b[0]) % self.q
        for i in range(1, self.n):
            b[i] = (-b[i]) % self.q

        qt  = self.q >> self.nut
        b = poly_rshift(b, self.nut, qt)

        #   --- 9.  return ( (vk := seed, b'), sk:= (vk, [[s]]) )
        vk = (seed, b)
        msk = (seed, b, ms_ntt)
        return msk, vk

    def sign_msg(self, msk, tr, msg):
        """Signing procedure of Plover (core: signs the mu hash)."""


        #   --- 1.  (vk, [[s]]) := [[sk]], (seed, b') := vk      [ caller ]
        (seed, b, ms_ntt) = msk

        #   --- 2.  mu := H( H(vk) || msg )                     [ caller ]

        #   --- 3.  A := ExpandA(seed)
        a_ntt = ntt(self._expand_a(seed))

        #   (restart position.)
        rsp_norms = False
        while not rsp_norms:

            #   --- x.  salt <- {0,1}^{2*bitsec}
            salt = self.random_bytes(self.salt_sz)

            u = self._msg_hash(salt, tr, msg)

            #   --- x. p <- 2 x ZeroEncoding()
            mp = [ self._vec_add_rep_noise( self.up, i, self.rep ) for i in range(2) ]

            #   --- x. [[w]] <- a * [[p_2]] + [[p_1]]
            mp_ntt = mat_ntt(mp)
            mw_ntt = [[ None ] for _ in range(self.d)]
            for i_d in range(self.d):
                mw_ntt[i_d] = poly_add(mp_ntt[0][i_d], 
                                       mul_ntt(a_ntt, mp_ntt[1][i_d]))
            self._refresh(mw_ntt)

            #   --- x. w <- Decode([[w]])
            w_ntt = self._decode(mw_ntt)
            w = intt(w_ntt)

            #   --- x. u' <- u - w
            u_ = poly_sub(u, w)

            #   --- x. u' -> q_1 * c_1 + c_2
            c1 = [0] * self.n

            # centers the distributions of c_1 and c_2
            mid1 = self.maxc1
            mid2 = 1 << (self.logdivide - 1)
            mid = mid1 * (1 << self.logdivide) + mid2

            for i in range(self.n):
                c1[i] = (((u_[i] + mid) % self.q) >> self.logdivide) - mid1

            c1_ntt = ntt(c1.copy())

            #   --- 12. [[s]] <- Refresh([[s]])
            self._refresh(ms_ntt)

            #   --- x.  [[z_2]] <- c_1 * [[s]] + [[p_2]]
            mz2_ntt = [ [ None ] for _ in range(self.d) ]
            for i_d in range(self.d):
                mz2_ntt[i_d] = poly_add(mul_ntt(c1_ntt, ms_ntt[i_d]),
                                        mp_ntt[1][i_d])
            self._refresh(mz2_ntt)
            z2_ntt = self._decode(mz2_ntt)
            z2 = intt(z2_ntt)

            #   --- 19. sig := (c_hash, h, z)                   [caller]

            #   --- 20. if CheckBounds(sig) = FAIL goto Line 4
            z2 = poly_center(z2)
            z3 = poly_center(c1)

            rsp_norms = self._check_bounds(u, a_ntt, b, z2, z3)

        #   --- 21. return sig
        sig = (salt, z2, z3)
        return sig


    def verify_msg(self, vk, tr, msg, sig):
        """Verification procedure of Plover (core: verifies msg)."""

        #   --- 1.  (c hash, h, z) := sig, (seed, b) := vk
        (salt, z2, z3) = sig
        (seed, b) = vk

        #   --- 2.  if CheckBounds(h, z) = FAIL return FAIL

        #   --- 4.  A := ExpandA(seed)
        a_ntt = ntt(self._expand_a(seed))

        #   --- 5. 
        u = self._msg_hash(salt, tr, msg)

        if self._check_bounds(u, a_ntt, b, z2, z3) == False:
            return False

        return True

    def set_random(self, random_bytes):
        """Set the key material RBG."""
        self.random_bytes   =   random_bytes

    def set_masking(self, masking_poly):
        """Set masking generator."""
        self.masking_poly = masking_poly

    #   --- internal methods ---

    def _compute_metrics(self):
        """Derive rejection bounds from parameters."""

        # max absolute value of a coefficient of c_1
        self.maxc1 = (((self.q-1) >> self.logdivide) + 1) >> 1

        sigp = sqrt(self.rep * ((1 << self.up)**2+1)/12)
        self.sq_beta = floor( 2**2 * (self.n * (2 * (sigp ** 2) + (self.maxc1 ** 2 + (1 << self.logdivide) ** 2) / 4) ) )
        # additional error introduced by the rounding of t
        self.sq_beta += self.n**2*((1<<(self.nut-1))*self.maxc1)**2

    def _check_bounds(self, u, a_ntt, t, z2, c1):
        """Check signature bounds. Return True iff bounds are acceptable."""

        #   this function only checks the norms; steps 1 and 2 are external.
        #   --- 1.  if |sig| != |sig|default return FAIL        [caller]
        #   --- 2.  (c hash, h, z) := sig                       [caller]

        t2 = poly_lshift(t, self.nut)
        t2 = [c - (1 << (self.nut-1)) for c in t2]
        t_ntt = ntt(t2.copy())

        #   --- 6.  z1 = u - a*z_2 + t*c_1
        z2_ntt = ntt(z2.copy())
        c1_ntt = ntt(c1.copy())

        z1_ntt = [a_ntt[i]*z2_ntt[i] + t_ntt[i]*c1_ntt[i] for i in range(self.n)]
        z1 = intt(z1_ntt)
        z1 = poly_sub(u, z1)

        # infinite norm for c1
        for c in c1:
            if abs(c) > self.maxc1:
                return False

        # norm 2 for the full vector
        sq_bound = self.sq_beta
        sq_n = 0
        for v in poly_center(z1):
            sq_n += v*v
        for v in poly_center(z2):
            sq_n += v*v
        for v in poly_center(c1):
            sq_n += v*v

        return sq_n <= sq_bound

    def _decode(self, mp):
        """Decode(): Collapse shares into a single polynomial."""
        r = mp[0].copy()
        for p in mp[1:]:
            r = poly_add(r, p)
        return r

    def _zero_encoding(self):
        """ZeroEncoding(): Create a masked encoding of zero."""

        z = [ [0] * self.n for _ in range(self.d) ]
        i = 1
        #   same ops as with recursion, but using nested loops
        while i < self.d:
            for j in range(0, self.d, 2 * i):
                for k in range(j, j + i):
                    r = self.masking_poly(self.n)
                    z[k] = poly_add(z[k], r)
                    z[k + i] = poly_sub(z[k + i], r)
            i <<= 1
        return z

    def _refresh(self, v):
        """Refresh(): Refresh shares via ZeroEncoding."""
        z = self._zero_encoding()
        for i in range(self.d):
            v[i] = poly_add(v[i], z[i])

    def _xof_sample_q(self, seed):
        """Expand a seed to n uniform values [0,q-1] using a XOF."""
        blen = (self.q_bits + 7) // 8
        mask = (1 << self.q_bits) - 1

        xof = SHAKE256.new(seed)
        v = [0] * self.n
        i = 0
        while i < self.n:
            z = xof.read(blen)
            x = int.from_bytes(z, BYTEORDER) & mask
            if (x < self.q):
                v[i] = x
                i += 1
        return v

    def _expand_a(self, seed):
        """ExpandA(): Expand "seed" into a polynomial."""

        #   XOF( 'A' || seed )
        xof_in  = bytes([ord('A'), 0, 0, 0, 0, 0, 0, 0]) + seed
        return self._xof_sample_q(xof_in)

    def _xof_sample_u(self, seed, u):
        """Sample a keyed uniform noise polynomial."""
        blen = (u + 7) // 8
        mask = (1 << u) - 1
        mid = (1 << u) // 2
        xof = SHAKE256.new(seed)
        r = [0] * self.n
        for i in range(self.n):
            z = xof.read(blen)
            x = int.from_bytes(z, BYTEORDER) & mask
            x ^= mid        # two's complement sign (1=neg)
            r[i] = (x - mid) % self.q
        return r


    def _vec_add_rep_noise(self, u, i_v, rep):
        """Repeatedly add uniform noise."""

        #   --- 1.  [[v]] <- (0_G)^d
        v = [[0] * self.n for j in range(self.d)]

        #   --- 2.  for i_rep in [rep] do
        for i_rep in range(rep):

            #   --- 3. for j in [d] do
            for j in range(self.d):

                #   --- 4.  rho <- {0,1}^lambda
                sigma = self.random_bytes(self.sec)

                #   --- 5.  hdr_u = ( 'u', rep, i, j, 0, 0, 0, 0 )
                hdr_u   = bytes([ord('u'), i_rep, i_v, j,
                                        0, 0, 0, 0]) + sigma

                #   --- 6.  v_i,j <- v_i,j + SampleU( hdr_u, sigma, u )
                r       = self._xof_sample_u(hdr_u, u)
                v[j] = poly_add(v[j], r)

            #   --- 7. Refresh([[v]])
            self._refresh(v)

        #   --- 8. Return [[v]]
        return v

    def _msg_hash(self, salt, tr, msg):
        """Compute the message hash for the signature (a single hash)."""

        xof_in  = bytes([ord('h'), 0, 0, 0, 0, 0, 0, 0]) + salt + tr + msg

        return self._xof_sample_q(xof_in)

#   --- some testing code ----------------------------------------------

if (__name__ == "__main__"):

    def chksum(v, q=549824583172097,g=15,s=31337):
        """Simple recursive poly/vector/matrix checksum routine."""
        if isinstance(v, int):
            return ((g * s + v) % q)
        elif isinstance(v, list):
            for x in v:
                s = chksum(x,q=q,g=g,s=s)
        return s

    def chkdim(v, s=''):
        t = v
        while isinstance(t, list):
            s += '[' + str(len(t)) + ']'
            t = t[0]
        s += ' = ' + str(chksum(v))
        return s

    #   one instance here for testing
    iut = PloverSign(  bitsec=128, q=PLOVERSIGN_Q, logdivide=37, nut=21, rep=8, ut=6,
                    up=42, n=PLOVERSIGN_N, d=1)

    #   initialize nist pseudo random
    entropy_input = bytes(range(48))
    drbg = NIST_KAT_DRBG(entropy_input)

    iut.set_random(drbg.random_bytes)
    iut.set_masking(MaskRandom().random_poly)

    print(f'name = {iut.name}')

    print("=== Keygen ===")
    msk, vk = iut.keygen()
    print(f"key: seed = {msk[0].hex().upper()}")
    print(chkdim(msk[1], 'key: t'))
    print(chkdim(msk[2], 'key: s'))

    for _ in range(5):
        print("=== Sign ===")
        tr = bytes(range(iut.tr_sz))
        msg = bytes(range(3))

        sig = iut.sign_msg(msk, tr, msg)
        print(f"sig: salt = {sig[0].hex().upper()}")
        print(chkdim(sig[1], 'sig: z2'))
        print(chkdim(sig[2], 'sig: c1'))

        print("=== Verify ===")
        rsp = iut.verify_msg(vk, tr, msg, sig)
        print(rsp)
        assert(rsp is True)

