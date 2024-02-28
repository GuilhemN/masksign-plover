//  plover_serial.c
//  Copyright (c) 2023 Plover Signature Team. See LICENSE.

//  === Plover signature scheme -- Serialize/deserialize.

#include <string.h>
#include "plover_core.h"
#include "plover_serial.h"
#include "plat_local.h"
#include "polyr.h"
#include "xof_sample.h"
#include "nist_random.h"
#include "mont64.h"
#include "sha3_t.h"

//  Encode vector v[PLOVER_N] as packed "bits" sized elements to  *b".
//  Return the number of bytes written -- at most ceil(PLOVER_N * bits/8).

static inline size_t inline_encode_bits(uint8_t *b, const int64_t v[PLOVER_N],
                                        size_t bits)
{
    size_t i, j, l;
    int64_t x, m;

    i = 0;  //  source word v[i]
    j = 0;  //  destination byte b[j]
    l = 0;  //  number of bits in x
    x = 0;  //  bit buffer

    m = (1llu << bits) - 1llu;

    while (i < PLOVER_N) {
        while (l < 8 && i < PLOVER_N) {
            x |= (v[i++] & m) << l;
            l += bits;
        }
        while (l >= 8) {
            b[j++] = (uint8_t)(x & 0xFF);
            x >>= 8;
            l -= 8;
        }
    }
    if (l > 0) {
        b[j++] = (uint8_t)(x & 0xFF);
    }
    return j;  //   return number of bytes written
}

//  Decode bytes from "*b" as PLOVER_N vector elements of "bits" each.
//  The decoding is unsigned if "is_signed"=false, two's complement
//  signed representation assumed if "is_signed"=true. Return the
//  number of bytes read -- upper bounded by ceil(PLOVER_N * bits/8).

static inline size_t inline_decode_bits(int64_t v[PLOVER_N], const uint8_t *b,
                                        size_t bits, bool is_signed)
{
    size_t i, j, l;
    int64_t x, m, s;

    i = 0;  //  source byte b[i]
    j = 0;  //  destination word v[j]
    l = 0;  //  number of bits in x
    x = 0;  //  bit buffer

    if (is_signed) {
        s = 1llu << (bits - 1);  // extract sign bit
        m = s - 1;
    } else {
        s = 0;  //  sign bit ignored
        m = (1llu << bits) - 1;
    }

    while (j < PLOVER_N) {

        while (l < bits) {
            x |= ((uint64_t)b[i++]) << l;
            l += 8;
        }
        while (l >= bits && j < PLOVER_N) {
            v[j++] = (x & m) - (x & s);
            x >>= bits;
            l -= bits;
        }
    }

    return i;  //   return number of bytes read
}

//  === Interface

//  Encode the public key "pk" to bytes "b". Return length in bytes.

size_t plover_encode_pk(uint8_t *b, const plover_pk_t *pk)
{
    size_t l;

    l = 0;  //  l holds the length

    //  encode A seed
    memcpy(b + l, pk->a_seed, PLOVER_AS_SZ);
    l += PLOVER_AS_SZ;

    //  encode t vector
    //  domain is q_t; has log2(q) - log(p_t) bits
    l += inline_encode_bits(b + l, pk->b, PLOVER_Q_BITS - PLOVER_NUT);
    return l;
}

//  Decode a public key from "b" to "pk". Return length in bytes.

size_t plover_decode_pk(plover_pk_t *pk, const uint8_t *b)
{
    size_t l;

    l = 0;

    //  decode A seed
    memcpy(pk->a_seed, b + l, PLOVER_AS_SZ);
    l += PLOVER_AS_SZ;

    //  decode t vector
    //  domain is q_t; has log2(q) - log(p_t) bits, unsigned
    l += inline_decode_bits(pk->b, b + l, PLOVER_Q_BITS - PLOVER_NUT, false);

    //  also set the tr field
    shake256(pk->tr, PLOVER_TR_SZ, b, l);

    return l;
}

//  Encode secret key "sk" to bytes "b". Return length in bytes.

size_t plover_encode_sk(uint8_t *b, const plover_sk_t *sk)
{
    size_t j, l;
    uint8_t buf[PLOVER_MK_SZ + 8];
    int64_t r[PLOVER_N], s0[PLOVER_N], e0[PLOVER_N];

    //  encode public key
    l = plover_encode_pk(b, &sk->pk);

    //  make a copy of share 0
    polyr_copy(s0, sk->s[0]);

    memset(buf, 0x00, 8);  //   domain header template
    buf[0] = 'K';

    //  shares 1, 2, ..., d-1
    for (j = 1; j < PLOVER_D; j++) {

        randombytes(b + l, PLOVER_MK_SZ);    //   key_j
        memcpy(buf + 8, b + l, PLOVER_MK_SZ);  // store in secret key
        l += PLOVER_MK_SZ;

        //  XOF( 'K' || share j || key_j )
        buf[1] = j;  // update domain header
        xof_sample_q(r, buf, PLOVER_MK_SZ + 8);
        polyr_subq(s0, s0, r);  //    s0 <- s0 - r
        polyr_addq(s0, s0, sk->s[j]);  //  s0 <- s0 + s_j
    }

    //  encode the zeroth share (in full)
    polyr_ntt_smul(s0, s0, MONT_R);
    l += inline_encode_bits(b + l, s0, PLOVER_Q_BITS);

    return l;
}

//  Decode secret key "sk" to bytes "b". Return length in bytes.

size_t plover_decode_sk(plover_sk_t *sk, const uint8_t *b)
{
    size_t j, l;
    uint8_t buf[PLOVER_MK_SZ + 8];

    //  decode public key
    l = plover_decode_pk(&sk->pk, b);

    memset(buf, 0x00, 8);  //   domain header template
    buf[0] = 'K';

    //  expand shares 1, 2, ..., d-1 from keys
    for (j = 1; j < PLOVER_D; j++) {
        //  copy key
        memcpy(buf + 8, b + l, PLOVER_MK_SZ);
        l += PLOVER_MK_SZ;

        //  XOF( 'K' || share j || key_j )
        buf[1] = j;  // update domain header
        xof_sample_q(sk->s[j], buf, PLOVER_MK_SZ + 8);
    }

    //  decode the zeroth share (in full)
    l += inline_decode_bits(sk->s[0], b + l, PLOVER_Q_BITS, false);

    return l;
}

//  macro for encoding n bits from y
//  (note -- returns from function on overflow)
#define ENC_SIG_PUT_BITS(y,n) { \
    while (n > 0) {             \
        n--;                    \
        z |= (y & 1) << k;      \
        y >>= 1;                \
        k++;                    \
        if (k == 8) {           \
            if (l >= b_sz)      \
                return 0;       \
            b[l++] = z;         \
            k = 0;              \
            z = 0;              \
        }                       \
    }                           \
}

/* see inner.h */

//  Encode signature "sig" to "*b" of max "b_sz" bytes. Return length in
//  bytes or zero in case of overflow.

size_t plover_encode_sig(uint8_t *buf, size_t b_sz, const plover_sig_t *sig)
{
    size_t u, v;
    int64_t acc;
    unsigned acc_len;

    acc = 0;
    acc_len = 0;
    v = 0;

    /**
     * Encode salt.
     */
    memcpy(buf, sig->salt, PLOVER_SALT_SZ);
    v += PLOVER_SALT_SZ;

    /**
     * First encode c_1.
     */
    int64_t c1[PLOVER_N];
    polyr_nonneg(c1, sig->c1, 1 << (PLOVER_SIG_BITS_C1 + 1)); // Complement 2
    v += inline_encode_bits(buf + v, sig->c1, PLOVER_SIG_BITS_C1 + 1);

    /*
     * We want to finish in a decent time, so the varying part of the encoding
     * of z_2 should not be too big.
     */
    int64_t maxval = (1L << PLOVER_SIG_RATE) * 16;
    for (u = 0; u < PLOVER_N; u++)
    {
        if (sig->z[u] < -maxval || sig->z[u] > maxval)
        {
            return 0;
        }
    }

    for (u = 0; u < PLOVER_N; u++)
    {
        int64_t t;
        uint64_t w;

        /*
         * Get sign and absolute value of next integer; push the
         * sign bit.
         */
        acc <<= 1;
        t = sig->z[u];
        if (t < 0)
        {
            t = -t;
            acc |= 1;
        }
        acc_len += 1; // push sign bit
        w = (uint64_t)t;

        /*
         * Push the low bits of the absolute value.
         */
        acc <<= PLOVER_SIG_RATE;
        acc |= w & ((1uL << PLOVER_SIG_RATE) - 1uL);
        w >>= PLOVER_SIG_RATE;

        acc_len += PLOVER_SIG_RATE;

        /*
         * Reduce accumulator.
         */
        while (acc_len >= 8)
        {
            acc_len -= 8;
            if (buf != NULL)
            {
                if (v >= b_sz)
                {
                    return 0;
                }
                buf[v] = (uint8_t)(acc >> acc_len);
            }
            v++;
        }

        /*
         * Push as many zeros as necessary, then a one. Since the
         * absolute value is bounded by an appropriate value, we can only range up
         * to 15 at this point, thus we will add at most 16 bits here which fits in 
         * the accumulator which previously contains at most 7 bits.
         */
        acc <<= (w + 1);
        acc |= 1;
        acc_len += w + 1;

        /*
         * Reduce accumulator.
         */
        while (acc_len >= 8)
        {
            acc_len -= 8;
            if (buf != NULL)
            {
                if (v >= b_sz)
                {
                    return 0;
                }
                buf[v] = (uint8_t)(acc >> acc_len);
            }
            v++;
        }
    }

    /*
     * Flush remaining bits (if any).
     */
    if (acc_len > 0)
    {
        if (buf != NULL)
        {
            if (v >= b_sz)
            {
                return 0;
            }
            buf[v] = (uint8_t)(acc << (8 - acc_len));
        }
        v++;
    }

    return v;
}

//  Encode signature "sig" to "*b" of max "b_sz" bytes. Return length in
//  bytes or zero in case of overflow.

#undef ENC_SIG_PUT_BITS

//  macro that gets a single bit
#define DEC_SIG_GET_BIT(bit) {  \
    bit = (z >> k) & 1;         \
    k++;                        \
    if (k == 8) {               \
        if (l >= b_sz)          \
            return 0;           \
        z = b[l++];             \
        k = 0;                  \
    }                           \
}

//  decode bytes "b" into signature "sig". Return length in bytes.

/* see inner.h */
// size_t Zf(comp_decode)(discrete_vector *x, const void *in, size_t max_in_len,
//                        unsigned rate)
size_t plover_decode_sig(plover_sig_t *sig, const uint8_t *buf, size_t b_sz)  
{
    size_t u, v;
    int64_t acc;
    unsigned acc_len;

    acc = 0;
    acc_len = 0;
    v = 0;

    /**
     * Decode salt.
     */
    memcpy(sig->salt, buf, PLOVER_SALT_SZ);
    v += PLOVER_SALT_SZ;

    // First decode c_1
    v += inline_decode_bits(sig->c1, buf + v, PLOVER_SIG_BITS_C1 + 1, true);

    uint64_t maxval =
        (1uL << PLOVER_SIG_RATE) * 16uL; // After removing the lowest rate bits, the value must
                          // not be too big (< 16), so that our computations fit
                          // in a 32 bits integer Also assumes that rate <= 8

    for (u = 0; u < PLOVER_N; u++)
    {
        uint64_t b, s, m;

        /*
         * Fill accumulator
         */
        while (acc_len < PLOVER_SIG_RATE + 1)
        {
            // not enough bits left in the accumulator
            if (v >= b_sz)
            {
                return 0;
            }

            acc = (acc << 8) | (uint32_t)buf[v++];
            acc_len += 8;
        }

        b = acc >> (acc_len - (PLOVER_SIG_RATE + 1));
        acc_len -= PLOVER_SIG_RATE + 1;
        s = b & (1uL << PLOVER_SIG_RATE);
        m = b & ((1uL << PLOVER_SIG_RATE) - 1uL);

        /*
         * Get next bits until a 1 is reached.
         */
        for (;;)
        {
            if (acc_len == 0)
            {
                if (v >= b_sz)
                {
                    return 0;
                }
                acc = (acc << 8) | (uint32_t)buf[v++];
                acc_len = 8;
            }
            acc_len--;
            if (((acc >> acc_len) & 1) != 0)
            {
                break;
            }
            m += 1L << PLOVER_SIG_RATE;
            if (m >= maxval)
            {
                return 0;
            }
        }

        /*
         * "-0" is forbidden.
         */
        if (s && m == 0)
        {
            return 0;
        }

        sig->z[u] = s ? -(int64_t)m : (int64_t)m;
    }

    /*
     * Unused bits in the last byte must be zero.
     */
    if ((acc & ((1uL << acc_len) - 1uL)) != 0)
    {
        return 0;
    }

    /**
     * Extra bytes of the signature must be zero.
     */
    for (; v < b_sz; v++) {
        if (buf[v] != 0) {
            return 0;
        }
    }

    return v;
}

#undef DEC_SIG_GET_BIT
