//  xof_sample.h
//  Copyright (c) 2023 Plover Signature Team. See LICENSE.

//  === Plover signature scheme -- Samplers and XOF functions

#ifndef _XOF_SAMPLE_H_
#define _XOF_SAMPLE_H_

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include "plover_param.h"

//  === Global namespace prefix
#ifdef PLOVER_
#define xof_sample_q    PLOVER_(xof_sample_q)
#define xof_sample_u    PLOVER_(xof_sample_u)
#define xof_chal_mu     PLOVER_(xof_chal_mu)
#define xof_msg_hash    PLOVER_(xof_msg_hash)
#define xof_chal_poly   PLOVER_(xof_chal_poly)
#endif

#ifdef __cplusplus
extern "C" {
#endif

//  Expand "seed" of "seed_sz" bytes to a uniform polynomial (mod q).
//  The input seed is assumed to already contain domain separation.
void xof_sample_q(int64_t r[PLOVER_N], const uint8_t *seed, size_t seed_sz);

//  Sample "bits"-wide signed coefficients from "seed[seed_sz]".
//  The input seed is assumed to alredy contain domain separation.
void xof_sample_u(int64_t r[PLOVER_N], int bits,
                  const uint8_t *seed, size_t seed_sz);

//  Hash salt with message to produce target polynomial "u".
void xof_msg_hash(int64_t cp[PLOVER_N], const uint8_t salt[PLOVER_SALT_SZ],
                  const uint8_t tr[PLOVER_TR_SZ], const uint8_t *m, size_t m_sz);

#ifdef __cplusplus
}
#endif

//  _XOF_SAMPLE_H_
#endif
