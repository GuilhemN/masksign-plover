//  plover_param.h
//  Copyright (c) 2023 Plover Signature Team. See LICENSE.

//  === Plover signature scheme -- Derived parameters.

#ifndef _PLOVER_PARAM_H_
#define _PLOVER_PARAM_H_

//  select a default parameter if somehow not defied
#if !defined(NIST_KAT) && !defined(BENCH_TIMEOUT)
#include "param_select.h"
#endif

//  include the parameter list
#include "param_list.h"

//  Byte size of symmetric keys / pre-image security
#define PLOVER_SEC    (PLOVER_KAPPA / 8)

//  Byte size for collision resistant hashes
#define PLOVER_CRH    ((2 * PLOVER_KAPPA) / 8)

//  Size of A_seed
#define PLOVER_AS_SZ  PLOVER_SEC

//  Size of public key hash used in BUFFing -- needs CRH
#define PLOVER_TR_SZ  PLOVER_CRH

//  Size of signture salt
#define PLOVER_SALT_SZ  PLOVER_CRH


//  Size of "mask keys" in serialized secret key
#define PLOVER_MK_SZ  PLOVER_SEC

//  shared / derived parameters
#if (PLOVER_Q == 2004477689857l) && (PLOVER_N == 2048)
#define PLOVER_Q_BITS 41
#define PLOVER_LGN    11
#else
#error  "No known parameter defined."
#endif

#define PLOVER_QMSK   ((1LL << PLOVER_Q_BITS) - 1)
#define PLOVER_QT     (PLOVER_Q >> PLOVER_NUT)

//  _PLOVER_PARAM_H_
#endif
