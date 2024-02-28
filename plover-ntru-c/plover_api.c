//  plover_api.c
//  Copyright (c) 2023 Plover Signature Team. See LICENSE.

//  === Plover signature scheme -- NIST KAT Generator API.

#include <string.h>

#include "api.h"
#include "plover_core.h"
#include "plover_serial.h"
#include "xof_sample.h"
#include "mont64.h"

//  Generates a keypair - pk is the public key and sk is the secret key.

int
crypto_sign_keypair(  unsigned char *pk, unsigned char *sk)
{
    plover_pk_t   r_pk;           //  internal-format public key
    plover_sk_t   r_sk;           //  internal-format secret key

    plover_core_keygen(&r_pk, &r_sk); //  generate keypair

    //  serialize
    if (CRYPTO_PUBLICKEYBYTES != plover_encode_pk(pk, &r_pk) ||
        CRYPTO_SECRETKEYBYTES != plover_encode_sk(sk, &r_sk))
        return -1;

    return 0;
}

//  Sign a message: sm is the signed message, m is the original message,
//  and sk is the secret key.

int
crypto_sign(unsigned char *sm, unsigned long long *smlen,
            const unsigned char *m, unsigned long long mlen,
            const unsigned char *sk)
{
    plover_sk_t   r_sk;           //  internal-format secret key
    plover_sig_t  r_sig;          //  internal-format signature
    size_t  sig_sz;
    //  deserialize secret key
    if (CRYPTO_SECRETKEYBYTES != plover_decode_sk(&r_sk, sk)) 
        return -1;

    //  several trials may be needed in case of signature size overflow
    do {
        plover_core_sign(&r_sig, &r_sk, m, mlen);          //  create signature

        //  The NIST API expects an "envelope" consisting of the message
        //  together with signature. we put the signature first.
        sig_sz = plover_encode_sig(sm, CRYPTO_BYTES, &r_sig);
    } while (sig_sz == 0);

    memset(sm + sig_sz, 0, CRYPTO_BYTES - sig_sz);  //  zero padding
    memcpy(sm + CRYPTO_BYTES, m, mlen);             //  add the message

    *smlen = mlen + CRYPTO_BYTES;

    return  0;
}

//  Verify a message signcoreature: m is the original message, sm is the signed
//  message, pk is the public key.

int
crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk)
{
    plover_pk_t   r_pk;           //  internal-format public key
    plover_sig_t  r_sig;          //  internal-format signature
    size_t      m_sz;

    //  deserialize public key, signature with a consistency check
    if (smlen < CRYPTO_BYTES ||
        CRYPTO_PUBLICKEYBYTES != plover_decode_pk(&r_pk, pk) ||
        CRYPTO_BYTES != plover_decode_sig(&r_sig, sm, CRYPTO_BYTES))
        return -1;
    m_sz = smlen - CRYPTO_BYTES;

    //  verification
    if (!plover_core_verify(&r_sig, &r_pk, sm + CRYPTO_BYTES, m_sz))
        return -1;

    //  store the length and move the "opened" message
    memcpy(m, sm + CRYPTO_BYTES, m_sz);
    *mlen = m_sz;

    return  0;
}

