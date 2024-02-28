//  api.h
//  === NIST Signature API

#ifndef _API_H_
#define _API_H_

#include "plover_param.h"

//  Set these three values apropriately for your algorithm
#define CRYPTO_SECRETKEYBYTES   PLOVER_SK_SZ
#define CRYPTO_PUBLICKEYBYTES   PLOVER_PK_SZ
#define CRYPTO_BYTES            PLOVER_SIG_SZ

// Change the algorithm name
#define CRYPTO_ALGNAME          PLOVER_NAME

int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk);

int
crypto_sign(unsigned char *sm, unsigned long long *smlen,
            const unsigned char *m, unsigned long long mlen,
            const unsigned char *sk);

int
crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk);

/* _API_H_ */
#endif
