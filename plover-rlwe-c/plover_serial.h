//  plover_serial.h
//  Copyright (c) 2023 Plover Signature Team. See LICENSE.

//  === Plover signature scheme -- Serialize/deserialize.

#ifndef _PLOVER_SERIAL_H_
#define _PLOVER_SERIAL_H_

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include "plover_param.h"

//  === Global namespace prefix

#ifdef PLOVER_
#define plover_encode_pk PLOVER_(encode_pk)
#define plover_decode_pk PLOVER_(decode_pk)
#define plover_encode_sk PLOVER_(encode_sk)
#define plover_decode_sk PLOVER_(decode_sk)
#define plover_encode_sig PLOVER_(encode_sig)
#define plover_decode_sig PLOVER_(decode_sig)
#endif

#ifdef __cplusplus
extern "C" {
#endif

//  Encode public key "pk" to bytes "b". Return length in bytes.
size_t plover_encode_pk(uint8_t *b, const plover_pk_t *pk);

//  Decode a public key from "b" to "pk". Return length in bytes.
size_t plover_decode_pk(plover_pk_t *pk, const uint8_t *b);

//  Encode secret key "sk" to bytes "b". Return length in bytes.
size_t plover_encode_sk(uint8_t *b, const plover_sk_t *sk);

//  Decode a secret key from "b" to "sk". Return length in bytes.
size_t plover_decode_sk(plover_sk_t *sk, const uint8_t *b);

//  Encode signature "sig" to "*b" of max "b_sz" bytes. Return length in
//  bytes or zero in case of overflow.
size_t plover_encode_sig(uint8_t *buf, size_t b_sz, const plover_sig_t *sig);

//  decode bytes "b" into signature "sig". Return length in bytes.
size_t plover_decode_sig(plover_sig_t *sig, const uint8_t *buf, size_t b_sz);

#ifdef __cplusplus
}
#endif

//  _PLOVER_SERIAL_H_
#endif
