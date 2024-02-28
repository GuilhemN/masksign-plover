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
                        q, logdivide, rep, ut, up, n, d,
                        masking_poly=MaskRandom().random_poly, 
                        masking_scalar=MaskRandom().uniform_q,
                        random_bytes=os.urandom, kappa=512):
        """Initialize a Plover instance."""

        self.name   =   f'Plover-{bitsec}-{d}'
        self.bitsec =   bitsec
        self.d      =   d
        self.q      =   q
        self.logdivide = logdivide
        self.q_bits =   self.q.bit_length()
        self.n      =   n
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
        self.masking_scalar = masking_scalar
        self.random_bytes = random_bytes

        #   calculate derived parmeters
        self._compute_metrics()

    def keygen(self):
        """Plover keypair generation."""

        #  --- 2.  [[f]] <- AddRepNoise([[f]], ut, rep)
        mf = self._vec_add_rep_noise( self.ut, 0, self.rep )
        mf_ntt = mat_ntt([mf])[0]

        #  --- 3.  while(e = 0)
        restart = True
        while restart:
            #  --- 4.   [[g]] <- AddRepNoise([[g]], ut, rep)
            mg = self._vec_add_rep_noise( self.ut, 0, self.rep )
            mg_ntt = mat_ntt([mg])[0]

            #  --- 5. [[g_\times]] <- PseudoInverse([[g]])
            mh_ntt = [[None]*self.n for _ in range(self.d)]

            restart = False
            for i in range(self.n):
                mfntt_i = [mf_ntt[i_d][i] for i_d in range(self.d)]
                mfntt_i[0] = (mfntt_i[0] - (1 << self.logdivide)) % self.q

                mgntt_i = [mg_ntt[i_d][i] for i_d in range(self.d)]

                # compute inverse of mgntt_i
                mgntt_i_inv = self.masked_pow(mgntt_i, PLOVERSIGN_PHIQ-1)

                # check if g[i] is invertible
                mg_check = self.masked_mul(mgntt_i, mgntt_i_inv)
                # decode g_check
                g_check = sum(mg_check) % self.q
                if g_check != 1:
                    restart = True
                    break
                
                #  --- 6.   [[h]] <- (divider-f)/g
                mhntt_i = self.masked_mul(mfntt_i, mgntt_i_inv)
                for i_d in range(self.d):
                    mh_ntt[i_d][i] = mhntt_i[i_d]
 
        h_ntt = self._decode(mh_ntt)
        h = intt(h_ntt)
        h = [(-hi) % self.q for hi in h]

        h_ntt = ntt(h.copy())

        #   --- 9.  return ( (vk := seed, h), sk:= (vk, [[g]]) )
        self._refresh(mg_ntt)

        vk = (h, )
        msk = (h, mg_ntt)

        return msk, vk

    def sign_msg(self, msk, tr, msg):
        """Signing procedure of Plover (core: signs the mu hash)."""

        (h, mg_ntt) = msk

        h_ntt = ntt(h.copy())

        #   (restart position.)
        rsp_norms = False
        while not rsp_norms:

            #   --- 1.  salt <- {0,1}^{2*kappa}
            salt = self.random_bytes(self.salt_sz)

            u = self._msg_hash(salt, tr, msg)

            #  --- 3.  [[p]] <- AddRepNoise([[p]], uw, rep)
            mp = [ self._vec_add_rep_noise( self.up, i, self.rep ) for i in range(2) ]

            #  --- 4. [[w]] <- a * [[p_2]] + [[p_1]]
            mp_ntt = mat_ntt(mp)
            mw_ntt = [[ None ] for _ in range(self.d)]
            for i_d in range(self.d):
                mw_ntt[i_d] = poly_add(mp_ntt[0][i_d], 
                                       mul_ntt(h_ntt, mp_ntt[1][i_d]))
            #   --- 5. w <- Unmask([[w]])
            w_ntt = self._decode(mw_ntt)
            w = intt(w_ntt)

            #   --- 6. c <- u - w
            u_ = poly_sub(u, w)

            #   --- 7. (c_1, c_2) <- Decompose(c)
            c1 = [0] * self.n

            # centers the distributions of c_1 and c_2
            mid1 = self.maxc1
            mid2 = 1 << (self.logdivide - 1)
            mid = mid1 * (1 << self.logdivide) + mid2

            for i in range(self.n):
                c1[i] = ((((u_[i] + mid) % self.q) >> self.logdivide) - mid1)%self.q

            c1_ntt = ntt(c1.copy())

            #   --- 8. [[g]] <- Refresh([[g]])
            self._refresh(mg_ntt)

            #   --- 9.  [[z_2]] <- c_1*[[g]] + [[p_2]]
            mz2_ntt = [ [ None ] for _ in range(self.d) ]
            for i_d in range(self.d):
                mz2_ntt[i_d] = poly_add(mul_ntt(c1_ntt, mg_ntt[i_d]),
                                        mp_ntt[1][i_d])
                
            #  --- 9.  z_2 := Unmask([[z_2]])
            z2_ntt = self._decode(mz2_ntt)

            # Compute z_1 from z_2, u and h
            h_ntt = ntt(h.copy())

            #  --- 10.  z_1' <- u-h*z_2
            z1_ = intt([(h_ntt[i]*z2_ntt[i]) % self.q for i in range(self.n) ])
            z1_ = poly_sub(u, z1_)

            z2 = intt(z2_ntt)
            z2 = poly_center(z2)

            #  --- 11.  if ||(z_1', z_2)|| > B, restart
            rsp_norms = self._check_bounds(z1_, z2)

        #   --- 21. return sig
        sig = (salt, z2)
        return sig


    def verify_msg(self, vk, tr, msg, sig):
        """Verification procedure of Plover (core: verifies mu)."""

        #   --- 1.  (c hash, h, z) := sig, (seed, t) := vk
        (salt, z2) = sig
        (h,) = vk

        h_ntt = ntt(h.copy())

        u = self._msg_hash(salt, tr, msg)

        #  --- 2.  z_1' <- u-h*z_2
        z2_ntt = ntt(z2.copy())
        z1_ = intt([(h_ntt[i]*z2_ntt[i]) % self.q for i in range(self.n) ])

        z1_ = poly_sub(u, z1_)

        #  --- 3.  if ||(z_1', z_2)|| > B, restart
        if self._check_bounds(z1_, z2) == False:
            return False

        return True

    def set_random(self, random_bytes):
        """Set the key material RBG."""
        self.random_bytes   =   random_bytes

    def set_masking(self, masking_poly):
        """Set masking generator."""
        self.masking_poly = masking_poly

    def set_masking_scalar(self, masking_scalar):
        """Set masking scalar generator."""
        self.masking_scalar = masking_scalar

    #   --- internal methods ---

    def _compute_metrics(self):
        """Derive rejection bounds from parameters."""

        # max absolute value of a coefficient of c_1
        self.maxc1 = (((self.q-1) >> self.logdivide) + 1) >> 1

        sigp = sqrt(self.rep * ((1 << self.up)**2+1)/12)
        self.sq_beta = floor( 2**2 * (self.n * (2 * (sigp ** 2) + ((1 << self.logdivide) ** 2) / 4) ) )

    def _check_bounds(self, z1, z2):
        """Check signature bounds. Return True iff bounds are acceptable."""

        #   this function only checks the norms; steps 1 and 2 are external.
        #   --- 1.  if |sig| != |sig|default return FAIL        [caller]
        #   --- 2.  (c hash, h, z) := sig                       [caller]

        # norm 2 for the full vector
        sq_bound = self.sq_beta
        sq_n = 0
        for v in poly_center(z1):
            sq_n += v*v
        for v in poly_center(z2):
            sq_n += v*v

        return sq_n <= sq_bound

    def _decode(self, mp):
        """Decode(): Collapse shares into a single polynomial."""
        self._refresh(mp)

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

    def masked_mul(self, a, b):
        """ 
        Implements masked multiplication.
        See https://eprint.iacr.org/2018/315.pdf
        """
        c = [(a[i] * b[i]) % self.q for i in range(self.d)]
        
        for i in range(self.d):
            for j in range(i+1, self.d):
                s = self.masking_scalar()
                s2 = (s + a[i]*b[j] + a[j]*b[i]) % self.q
                
                c[i] = (c[i] - s) % self.q
                c[j] = (c[j] + s2) % self.q

        return c
    
    def masked_pow(self, x, e):
        """ Computes x**e (mod q) when x is masked. """
        
        y = [1] + [0]*(self.d-1)
        while e > 0:
            if e & 1 == 1:
                y = self.masked_mul(y, x)
            x = self.masked_mul(x, x)
            e >>= 1

        return y

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
    iut = PloverSign(  bitsec=128, q=PLOVERSIGN_Q, logdivide=41, rep=8, ut=6,
                    up=42, n=PLOVERSIGN_N, d=4)

    #   initialize nist pseudo random
    entropy_input = bytes(range(48))
    drbg = NIST_KAT_DRBG(entropy_input)

    iut.set_random(drbg.random_bytes)
    iut.set_masking(MaskRandom().random_poly)
    iut.set_masking_scalar(MaskRandom().uniform_q)

    print(f'name = {iut.name}')

    print("=== Keygen ===")
    msk, vk = iut.keygen()
    print(f"key: ")
    print(chkdim(msk[0], 'key: f'))
    print(chkdim(msk[1], 'key: g'))

    for _ in range(5):
        print("=== Sign ===")
        tr = bytes(range(iut.tr_sz))
        msg = bytes(range(3))

        sig = iut.sign_msg(msk, tr, msg)
        print(f"sig: salt = {sig[0].hex().upper()}")
        print(chkdim(sig[1], 'sig: z'))

        print("=== Verify ===")
        rsp = iut.verify_msg(vk, tr, msg, sig)
        print(rsp)
        assert(rsp is True)

