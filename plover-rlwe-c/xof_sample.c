//  xof_sample.c
//  Copyright (c) 2023 Plover Signature Team. See LICENSE.

//  === Plover signature scheme -- Samplers and XOF functions

#include <string.h>

#include "plover_param.h"
#include "xof_sample.h"
#include "sha3_t.h"
#include "mont64.h"

//  Expand "seed" of "seed_sz" bytes to a uniform polynomial (mod q).
//  The input seed is assumed to already contain domain separation.

void xof_sample_q(int64_t r[PLOVER_N], const uint8_t *seed, size_t seed_sz)
{
    size_t i;
    int64_t x;
    uint8_t buf[8];
    sha3_t kec;

    //  sample from squeezed output
    sha3_init(&kec, SHAKE256_RATE);
    sha3_absorb(&kec, seed, seed_sz);
    sha3_pad(&kec, SHAKE_PAD);

    memset(buf, 0, sizeof(buf));
    for (i = 0; i < PLOVER_N; i++) {
        do {
            sha3_squeeze(&kec, buf, (PLOVER_Q_BITS + 7) / 8);
            x = get64u_le(buf) & PLOVER_QMSK;
        } while (x >= PLOVER_Q);
        r[i] = x;
    }
}

//  Sample "bits"-wide signed coefficients from "seed[seed_sz]".
//  The input seed is assumed to alredy contain domain separation.

void xof_sample_u(int64_t r[PLOVER_N], int bits,
                  const uint8_t *seed, size_t seed_sz)
{
    size_t i, blen;
    int64_t x, mask, mid;
    uint8_t buf[8];
    sha3_t kec;

    blen = (bits + 7) / 8;
    mask = (1ll << bits) - 1;
    mid = 1ll << (bits - 1);

    //  absorb seed
    sha3_init(&kec, SHAKE256_RATE);
    sha3_absorb(&kec, seed, seed_sz);
    sha3_pad(&kec, SHAKE_PAD);

    //  sample from squeezed outpu
    memset(buf, 0, sizeof(buf));
    for (i = 0; i < PLOVER_N; i++) {
        sha3_squeeze(&kec, buf, blen);
        x = get64u_le(buf) & mask;
        x ^= mid;  //   two's complement sign bit: 0=pos, 1=neg
        r[i] = mont64_cadd(x - mid, PLOVER_Q);
    }
}

//  Hash "w" vector with "mu" to produce challenge hash "ch".

void xof_msg_hash(int64_t cp[PLOVER_N], const uint8_t salt[PLOVER_SALT_SZ],
                  const uint8_t tr[PLOVER_TR_SZ], const uint8_t *m, size_t m_sz)
{
    sha3_t kec;
    uint8_t head[8], buf[8];
    size_t i;
    int64_t x;

    sha3_init(&kec, SHAKE256_RATE);         //  header
    head[0] = 'h';
    memset(head + 1, 0x00, 7);
    sha3_absorb(&kec, head, 8);

    sha3_absorb(&kec, salt, PLOVER_SALT_SZ); //  add salt
    sha3_absorb(&kec, tr, PLOVER_TR_SZ);     //  add hash of public key
    sha3_absorb(&kec, m, m_sz);            //  add message

    sha3_pad(&kec, SHAKE_PAD);

    memset(buf, 0, sizeof(buf));
    for (i = 0; i < PLOVER_N; i++)
    {
        do
        {
            sha3_squeeze(&kec, buf, (PLOVER_Q_BITS + 7) / 8);
            x = get64u_le(buf) & PLOVER_QMSK;
        } while (x >= PLOVER_Q);
        cp[i] = x;
    }

    sha3_clear(&kec);
}