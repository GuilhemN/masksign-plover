"""
plover_api.py
Copyright (c) 2023 Plover Signature Team. See LICENSE.

=== Masked Plover signature scheme: Serialization, parameters, BUFF interface.
"""

from Crypto.Hash import SHAKE256
from nist_kat_drbg import NIST_KAT_DRBG
from mask_random import MaskRandom
from ploversign_core import PloverSign
from polyr import *
from math import log, ceil

#   Encoding and decoding methods for NIST Test Vectors

class NIST_PloverSign(PloverSign):

    def __init__(self, sig_rate, sig_sz, *args, **kwargs):
        """This is a subclass that provides serialization and BUFF."""
        super().__init__(*args, **kwargs)

        #   nist serialization sizes
        self.pk_sz  =   (self.n * self.q_bits) // 8
        self.sk_sz  =   (self.pk_sz + (self.d - 1) * self.mk_sz +
                            (self.n * self.q_bits) // 8)

        self.sig_sz = sig_sz
        self.sig_rate = sig_rate

    @staticmethod
    def _encode_bits(v, bits):
        """Encode vector v of integers into bytes, 'bits' per element."""
        x = 0                           # bit buffer
        l = 0                           # number of bits in x
        i = 0                           # index in vector v[i]
        b = b''                         # zero-length array of bytes
        m = (1 << bits) - 1             # bit mask

        while i < len(v):
            while l < 8 and i < len(v):
                x |= (v[i] & m) << l    # load an integer into x
                i += 1
                l += bits
            while l >= 8:
                b += bytes([x & 0xFF])  # store a bytes from x
                x >>= 8
                l -= 8
        if l > 0:
            b += bytes([x])             # a byte with leftover bits

        return b

    """
    #   this is functionally equivalent but slower -- O(n^2)!
    @staticmethod
    def _encode_bits(v, bits):
        x = 0                   # bit buffer
        m = (1 << bits) - 1     # bit mask; "bits" ones
        for i in range(len(v)):
            x |= (v[i] & m) << (bits * i)
        return x.to_bytes( (bits * len(v) + 7) // 8, byteorder='little' )
    """

    @staticmethod
    def _decode_bits(b, bits, n, is_signed=False):
        """
        Decode bytes from 'b' into a vector of 'n' integers, 'bits' each.
        """
        x = 0                           # bit buffer
        i = 0                           # source byte index b[i]
        v = []                          # zero-length result vector
        l = 0                           # number of bits in x

        if is_signed:
            s = 1 << (bits - 1)         # sign bit is negative
            m = s - 1                   # mask bits-1 bits
        else:
            s = 0                       # unsigned: no sign bit
            m = (1 << bits) - 1         # mask given number of bits

        while len(v) < n:
            while l < bits:             # read bytes until full integer
                x |= int(b[i]) << l
                i += 1
                l += 8
            while l >= bits and len(v) < n: # write integer(s)
                v += [ (x & m) - (x & s) ]
                x >>= bits
                l -= bits

        return v, i     #   return the vector and number of bytes read

    def encode_pk(self, vk):
        """Serialize the signature verification (public) key."""
        (h,) = vk
        return self._encode_bits(h, self.q_bits)

    def encode_sk(self, msk):
        """Serialize the masked signing key."""
        (h, mg_ntt) = msk

        #   encode public key
        b = self.encode_pk((h,))

        #   copy share 0
        g0 = mg_ntt[0].copy()

        #   encode keys for shares 1, 2, ..., d-1
        for j in range(1, self.d):

            #   key_j for share j
            key =   self.random_bytes(self.mk_sz)
            b   +=  key

            #   update share 0
            #   XOF( 'K' || index j || (0 pad) || key_j )
            xof_in  = bytes([ord('K'), j, 0, 0, 0, 0, 0, 0]) + key

            r    = self._xof_sample_q(xof_in)
            g0   = poly_sub(g0, r)
            g0   = poly_add(g0, mg_ntt[j])

        #   encode share 0
        b += self._encode_bits(g0, self.q_bits)

        return b

    def decode_pk(self, b):
        """Decode the verification key from bytes."""
        h,pl = self._decode_bits(b, self.q_bits, self.n);
        l = pl
        vk = (h,)

        #   compute the "tr" hash from serialized public key
        tr = SHAKE256.new(b[0:l]).read(self.tr_sz)

        return vk, tr, l

    def decode_sk(self, b):
        """Decode a signing key from bytes."""

        #   decode public key
        vk, tr, l = self.decode_pk(b)
        h, = vk

        #   expand shares 1, 2, .. d-1
        mg_ntt = [None for _ in range(self.d)]

        for j in range(1, self.d):
            key =   b[l:l+self.mk_sz]
            l   +=  self.mk_sz

            #   XOF( 'K' || index j || (0 pad)  || key_j )
            xof_in  = bytes([ord('K'), j, 0, 0, 0, 0, 0, 0]) + key
            mg_ntt[j]    = self._xof_sample_q(xof_in)

        #   decode share zero
        mg_ntt[0],sl =   self._decode_bits(b[l:], self.q_bits, self.n)
        l           +=  sl

        msk = (h, mg_ntt)
        return msk, tr, l

    def _compress(self, v, slen):
        """
        Take as input a list of integers v and a bytelength slen, and
        return a bytestring of length slen that encode/compress v.
        If this is not possible, return False.

        For each coefficient of v:
        - the sign is encoded on 1 bit
        - the 7 lower bits are encoded naively (binary)
        - the high bits are encoded in unary encoding
        """
        u = ""
        for coef in v:
            # Encode the sign
            s = "1" if coef < 0 else "0"
            # Encode the low bits
            s += format((abs(coef) % (1 << self.sig_rate)), f'#0{self.sig_rate+2}b')[2:]
            # Encode the high bits
            s += "0" * (abs(coef) >> self.sig_rate) + "1"
            u += s
        # The encoding is too long
        if len(u) > 8 * slen:
            return False
        u += "0" * (8 * slen - len(u))
        w = [int(u[8 * i: 8 * i + 8], 2) for i in range(len(u) // 8)]
        x = bytes(w)
        return x


    def _decompress(self, x, slen, n):
        """
        Take as input an encoding x, a bytelength slen and a length n, and
        return a list of integers v of length n such that x encode v.
        If such a list does not exist, the encoding is invalid and we output False.
        """
        if (len(x) > slen):
            print("Too long")
            return False
        w = list(x)
        u = ""
        for elt in w:
            u += bin((1 << 8) ^ elt)[3:]
        v = []

        # Remove the last bits
        while u[-1] == "0":
            u = u[:-1]

        try:
            while (u != "") and (len(v) < n):
                # Recover the sign of coef
                sign = -1 if u[0] == "1" else 1
                # Recover the first low bits of abs(coef)
                low = int(u[1:1+self.sig_rate], 2)
                i, high = self.sig_rate+1, 0
                # Recover the high bits of abs(coef)
                while (u[i] == "0"):
                    i += 1
                    high += 1
                # Compute coef
                coef = sign * (low + (high << self.sig_rate))
                # Enforce a unique encoding for coef = 0
                if (coef == 0) and (sign == -1):
                    return False
                # Store intermediate results
                v += [coef]
                u = u[i + 1:]
            # In this case, the encoding is invalid
            if (len(v) != n):
                return False
            return v
        # IndexError is raised if indices are read outside the table bounds
        except IndexError:
            return False

    def encode_sig(self, sig):
        """Serialize a signature as bytes. No zero padding / length check."""
        (salt, z2) = sig
        s = salt                        #   bit string
        s += self._compress(z2, self.sig_sz - len(s))

        return s

    def decode_sig(self, s):
        """Deserialize a signature."""
        salt = s[0:self.salt_sz]
        i = self.salt_sz

        z2 = self._decompress(s[i:], self.sig_sz - i, self.n)

        return (salt, z2)

    #   interface that directly uses byte sequences

    def byte_keygen(self):
        """(API) Key pair generation directly into bytes."""
        msk, vk = self.keygen()
        return self.encode_pk(vk), self.encode_sk(msk)

    def byte_signature(self, msg, sk):
        """Detached signature generation directly from/to bytes."""
        msk, tr, _ = self.decode_sk(sk)
        sig_b = []
        while len(sig_b) != self.sig_sz:
            sig = self.sign_msg(msk, tr, msg)
            sig_b = self.encode_sig(sig)
            if len(sig_b) < self.sig_sz:
                sig_b += bytes([0] * (self.sig_sz - len(sig_b)))
        return sig_b

    def byte_verify(self, msg, sm, pk):
        """Detached Signature verification directly from bytes."""
        if len(sm) < self.sig_sz:
            return False
        vk, tr, _ = self.decode_pk(pk)
        sig = self.decode_sig(sm[0:self.sig_sz])
        return self.verify_msg(vk, tr, msg, sig)

    def byte_sign(self, msg, sk):
        """(API) Signature "envelope" generation directly from/to bytes."""
        sig = self.byte_signature(msg, sk)
        return sig + msg

    def byte_open(self, sm, pk):
        """(API) Signature verification directly from bytes."""
        msg = sm[self.sig_sz:]
        return self.byte_verify(msg, sm, pk), msg

### Instantiate Parameter sets

############################
### 128 bits of security ###
############################

# generated with params/gen_concrete.py

plover_128_1  = NIST_PloverSign(  bitsec=128, q=PLOVERSIGN_Q, logdivide=37, rep=8, ut=27, 
				up=36, n=PLOVERSIGN_N, d=1, sig_rate=37, sig_sz=14000)

plover_128_2  = NIST_PloverSign(  bitsec=128, q=PLOVERSIGN_Q, logdivide=37, rep=4, ut=27, 
				up=36, n=PLOVERSIGN_N, d=2, sig_rate=37, sig_sz=14000)

plover_128_4  = NIST_PloverSign(  bitsec=128, q=PLOVERSIGN_Q, logdivide=37, rep=2, ut=27, 
				up=36, n=PLOVERSIGN_N, d=4, sig_rate=37, sig_sz=14000)

plover_128_8  = NIST_PloverSign(  bitsec=128, q=PLOVERSIGN_Q, logdivide=37, rep=4, ut=26, 
				up=35, n=PLOVERSIGN_N, d=8, sig_rate=37, sig_sz=14000)

plover_128_16  = NIST_PloverSign(  bitsec=128, q=PLOVERSIGN_Q, logdivide=37, rep=2, ut=26, 
				up=35, n=PLOVERSIGN_N, d=16, sig_rate=37, sig_sz=14000)

plover_128_32  = NIST_PloverSign(  bitsec=128, q=PLOVERSIGN_Q, logdivide=37, rep=4, ut=25, 
				up=34, n=PLOVERSIGN_N, d=32, sig_rate=37, sig_sz=14000)

# end generated

plover_all = [
    plover_128_1, plover_128_2, plover_128_4, plover_128_8,
    plover_128_16, plover_128_32
]

p_test = plover_128_1

msk, pk = p_test.keygen()

# check pk encoding
pk2, tr, _ = p_test.decode_pk(p_test.encode_pk(pk))
assert(pk == pk2)

# check sk encoding
msk2, _, _ = p_test.decode_sk(p_test.encode_sk(msk))
assert(msk2[0] == msk[0]) # eq of h
assert(p_test._decode(msk2[1]) == p_test._decode(msk[1])) # eq g_ntt

# check signature encoding
msg = bytes(range(3))
sig = p_test.sign_msg(msk, tr, msg)
sig2 = p_test.decode_sig(p_test.encode_sig(sig))
assert(sig[0] == sig2[0]) # same salt
assert(sig[1] == sig2[1]) # same z

assert(p_test.verify_msg(pk, tr, msg, sig))

print("pk_sz", p_test.pk_sz)
print("sk_sz", p_test.sk_sz)