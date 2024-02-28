//  mont64.h
//  Copyright (c) 2023 Plover Signature Team. See LICENSE.

//  === Portable 64-bit Montgomery arithmetic

#ifndef _MONT64_H_
#define _MONT64_H_

#include "plat_local.h"
#include "plover_param.h"

// file generated with params/concrete-param-ring.py

#if (PLOVER_N != 2048 || PLOVER_Q != 2004477689857l)
#error "Unrecognized polynomial parameters N, Q"
#endif

/*
    n   = 2048
    q1  = 2**20 - 9*2**13
    q2  = 2**21 - 5*2**13
    q   = q1*q2
    r   = 2^64 % q
    rr  = r^2 % q
    ni  = lift(rr * Mod(n,q)^-1)
    qi  = lift(Mod(-q,2^64)^-1)
*/

//  Montgomery constants. These depend on Q and N
#define MONT_R 932779627440L
#define MONT_RR 1520476760420L
#define MONT_NI 1656785511718L
#define MONT_QI 2306286626255290367L        

// end generated      


//  Addition and subtraction

static inline int64_t mont64_add(int64_t x, int64_t y)
{
    return x + y;
}

static inline int64_t mont64_sub(int64_t x, int64_t y)
{
    return x - y;
}
//  Conditionally add m if x is negative

static inline int64_t mont64_cadd(int64_t x, int64_t m)
{
    int64_t t, r;

    XASSUME(x >= -m && x < m);

    t = x >> 63;
    r = x + (t & m);

    XASSERT(r >= 0 && r < m);
    XASSERT(r == x || r == x + m);

    return r;
}

//  Conditionally subtract m if x >= m

static inline int64_t mont64_csub(int64_t x, int64_t m)
{
    int64_t t, r;

    XASSUME(x >= 0 && x < 2 * m);
    XASSUME(m > 0);

    t = x - m;
    r = t + ((t >> 63) & m);

    XASSERT(r >= 0 && r < m);
    XASSERT(r == x || r == x - m);

    return r;
}

//  Montgomery reduction. Returns r in [-q,q-1] so that r == (x/2^64) mod q.

static inline int64_t mont64_redc(__int128 x)
{
    int64_t r;

    //  prove these input bounds
    XASSUME(x >= -(((__int128)1) << 111));
    XASSUME(x < (((__int128)1) << 111));

    r = x * MONT_QI;
    r = (x + ((__int128)r) * ((__int128)PLOVER_Q)) >> 64;

    //  prove that only one coditional addition is required
    XASSERT(r >= -PLOVER_Q && r < PLOVER_Q);

#ifdef XDEBUG
    //  this modular reduction correctness proof is too slow for SAT
    XASSERT(((((__int128)x) - (((__int128)r) << 64)) %
            ((__int128_t)PLOVER_Q)) == 0);
#endif
    return r;
}

//  Montgomery multiplication. r in [-q,q-1] so that r == (a*b)/2^64) mod q.

static inline int64_t mont64_mulq(int64_t x, int64_t y)
{
    int64_t r;

    r = mont64_redc(((__int128)x) * ((__int128)y));

    return r;
}

//  same with addition

static inline int64_t mont64_mulqa(int64_t x, int64_t y, int64_t z)
{
    int64_t r;

    r = mont64_redc(((__int128)x) * ((__int128)y) + ((__int128)z));

    return r;
}

//  _MONT64_H_
#endif
