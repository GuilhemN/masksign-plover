//  plover_core.h
//  Copyright (c) 2023 Plover Signature Team. See LICENSE.

//  === Plover signature scheme -- Core internal API.

#ifndef _PLOVER_CORE_H_
#define _PLOVER_CORE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include "plover_param.h"

//  === Global namespace prefix
#ifdef PLOVER_
#define plover_core_keygen PLOVER_(core_keygen)
#define plover_core_sign PLOVER_(core_sign)
#define plover_core_verify PLOVER_(core_verify)
#endif

//  === Internal structures ===

//  Plover public key
typedef struct {
    int64_t h[PLOVER_N];                      //  public key
    uint8_t tr[PLOVER_TR_SZ];                 //  hash of serialized public key
} plover_pk_t;

//  Plover secret key
typedef struct {
    plover_pk_t pk;                           //  copy of public key
    int64_t g_ntt[PLOVER_D][PLOVER_N];          //  d-masked secret key
} plover_sk_t;

//  Plover signature
typedef struct {
    uint8_t salt[PLOVER_SALT_SZ];             //  Signature salt
    int64_t z2[PLOVER_N];                      //  signature data
} plover_sig_t;

//  === Core API ===

//  Generate a public-secret keypair ("pk", "sk").
void plover_core_keygen(plover_pk_t *pk, plover_sk_t *sk);

//  Create a detached signature "sig" for digest "mu" using secret key "sk".
void plover_core_sign(plover_sig_t *sig, plover_sk_t *sk, const uint8_t *m,
                    size_t m_sz);

//  Verify that the signature "sig" is valid for digest "mu".
//  Returns true iff signature is valid, false if not valid.
bool plover_core_verify(const plover_sig_t *sig, const plover_pk_t *pk, const uint8_t *m,
                      size_t m_sz);

#ifdef __cplusplus
}
#endif

//  _PLOVER_CORE_H_
#endif
