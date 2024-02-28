//  plover_core.c
//  Copyright (c) 2023 Plover Signature Team. See LICENSE.

//  === Plover signature scheme -- core scheme.

#include <string.h>

#include "plat_local.h"
#include "plover_core.h"
#include "polyr.h"
#include "mont64.h"
#include "ct_util.h"
#include "xof_sample.h"
#include "nist_random.h"
#include "mask_random.h"


//  ZeroEncoding(d) -> [[z]]d
//  in-place version

static void zero_encoding(int64_t z[PLOVER_D][PLOVER_N], mask_random_t *mrg)
{
#if PLOVER_D == 1
    (void) mrg;
    polyr_zero(z[0]);
#else
    int i, j, d;
    int64_t r[PLOVER_N];

    //  d = 2
    for (i = 0; i < PLOVER_D; i += 2) {
        mask_random_poly(mrg, z[i], i);
        polyr_negm(z[i + 1], z[i], PLOVER_Q);
    }

    //  d = 4, 8, ..
    d = 2;
    while (d < PLOVER_D) {
        for (i = 0; i < PLOVER_D; i += 2 * d) {
            for (j = i; j < i + d; j++) {
                mask_random_poly(mrg, r, j);
                polyr_addq(z[j], z[j], r);
                polyr_subq(z[j + d], z[j + d], r);
            }
        }
        d <<= 1;
    }
#endif
}

//  Refresh([[x]]) -> [[x]]′

static void plover_refresh(int64_t x[PLOVER_D][PLOVER_N], mask_random_t *mrg)
{
#if PLOVER_D == 1
    (void) x;
    (void) mrg;
#else
    int i;
    int64_t z[PLOVER_D][PLOVER_N];

    //  --- 1.  [[z]] <- ZeroEncoding(d)
    zero_encoding(z, mrg);

    //  --- 2.  return [[x]]' := [[x]] + [[z]]
    for (i = 0; i < PLOVER_D; i++) {
        polyr_addq(x[i], x[i], z[i]);
    }
#endif
}

//  Refresh([[x]]) -> [[x]]′ ( NTT domain )

static void plover_ntt_refresh(int64_t x[PLOVER_D][PLOVER_N], mask_random_t *mrg)
{
#if PLOVER_D == 1
    (void) x;
    (void) mrg;
#else
    int i;
    int64_t z[PLOVER_D][PLOVER_N];

    //  --- 1.  [[z]] <- ZeroEncoding(d)
    zero_encoding(z, mrg);

    //  --- 2.  return [[x]]' := [[x]] + [[z]]
    for (i = 0; i < PLOVER_D; i++) {
#ifdef POLYR_Q32
        polyr2_split(z[i]);
#endif
        polyr_ntt_addq(x[i], x[i], z[i]);
    }
#endif
}

//  ExpandA(): Use domain separated XOF to create matrix elements

static void expand_a( int64_t a[PLOVER_N], const uint8_t seed[PLOVER_AS_SZ])
{
    uint8_t buf[PLOVER_AS_SZ + 8];

    //  --- 3.  hdrA := Ser8(65, 0, 0, 0, 0, 0, 0, 0)
    buf[0] = 'A';       //  ascii 65
    memset(buf + 1, 0x00, 8 - 1);

    //  --- 4.  Ai,j <- SampleQ(hdrA, seed)
    memcpy(buf + 8, seed, PLOVER_AS_SZ);
    xof_sample_q(a, buf, PLOVER_AS_SZ + 8);

    //  converted to NTT domain
    polyr_fntt(a);
}

//  Decode(): Collapse shares

static void plover_decode(int64_t r[PLOVER_N], const int64_t m[PLOVER_D][PLOVER_N], mask_random_t *mrg)
{
    plover_refresh(m, mrg);

#if PLOVER_D == 1
    polyr_copy(r, m[0]);
#else
    int i;

    polyr_addq(r, m[0], m[1]);
    for (i = 2; i < PLOVER_D; i++) {
        polyr_addq(r, r, m[i]);
    }
#endif
}

//  Decode(): Collapse shares (possibly split CRT arithmetic)

static void plover_ntt_decode(int64_t r[PLOVER_N], const int64_t m[PLOVER_D][PLOVER_N], mask_random_t *mrg)
{
    plover_refresh(m, mrg);

#if PLOVER_D == 1
    polyr_copy(r, m[0]);
#else
    int i;

    polyr_ntt_addq(r, m[0], m[1]);
    for (i = 2; i < PLOVER_D; i++) {
        polyr_ntt_addq(r, r, m[i]);
    }
#endif
}

//  AddRepNoise([[v]], u, rep) -> [[v]]
//  Add repeated noise to a polynomial (vector at index i_v)

static void add_rep_noise(  int64_t vi[PLOVER_D][PLOVER_N],
                            int i_v, int u, mask_random_t *mrg)
{
    int i_rep, i, j;
    uint8_t buf[PLOVER_SEC + 8];
    int64_t r[PLOVER_N];

    //  --- 1.  [[v]] <- (0_G)^d
    for (j = 0; j < PLOVER_D; j++) {
        for (i = 0; i < PLOVER_N; i++) {
            vi[j][i] = 0;
        }
    }

    //  --- 2.  for i_rep in [rep] do
    for (i_rep = 0; i_rep < PLOVER_REP; i_rep++) {

        //  --- 3.  for j in [d] do:
        for (j = 0; j < PLOVER_D; j++) {

            //  --- 4.  sigma <- {0,1}^kappa
            randombytes(buf + 8, PLOVER_SEC);

            //  --- 5.  hdr_u := Ser8('u' || i_rep || i_v || j || (0) || seed)
            buf[0] = 'u';       //  ascii 117
            buf[1] = i_rep;
            buf[2] = i_v;
            buf[3] = j;
            memset(buf + 4, 0x00, 8 - 4);

            //  --- 6.  v_ij <- v_ij + SampleU(hdr_u, sigma, u)
            xof_sample_u(r, u, buf, PLOVER_SEC + 8);
            polyr_addq(vi[j], vi[j], r);
        }

        //  --- [[v_i]] <- Refresh([[v_i]])
        plover_refresh(vi, mrg);
    }
}

//  "rounding" shift right

static inline void round_shift_r(int64_t *r, int64_t q, int s)
{
    int i;
    int64_t x, rc;

    rc = 1ll << (s - 1);
    for (i = 0; i < PLOVER_N; i++) {
        x = (r[i] + rc) >> s;
        r[i] = mont64_csub(x, q);
    }
}

//  CheckBounds(sig) -> {OK or FAIL}

static bool plover_check_bounds( const int64_t u[PLOVER_N], const int64_t a[PLOVER_N],
    const int64_t b[PLOVER_N], const plover_sig_t *sig)
{
    int i;
    int64_t z1[PLOVER_N], z2[PLOVER_N], c1[PLOVER_N], t[PLOVER_N];

    // Check infinity norm on c_1
    int64_t maxc1 = (((PLOVER_Q - 1) >> PLOVER_LOGD) + 1) >> 1;
    for (i = 0; i < PLOVER_N; i++) {
        if (sig->c1[i] > maxc1 || sig->c1[i] < -maxc1) {
            return false;
        }
    }

    // Recompute z_1
    polyr_copy(z2, sig->z);
    polyr_fntt(z2);

    polyr_copy(c1, sig->c1);
    polyr_fntt(c1);

    polyr_shlm(t, b, PLOVER_NUT, PLOVER_Q); //  .. - p_t * t ..
    polyr_fntt(t);

    // add a*z_2+b*c_1
    polyr_ntt_cmul(z1, a, z2);
    polyr_ntt_mula(z1, t, c1, z1);

    // z1 <- u - a*z_2 + b*c_1
    polyr_intt(z1);

    polyr_subq(z1, u, z1);

    // Check norm on vector z
    polyr_center(z1, z1, PLOVER_Q);
    __int128_t n = 0;
    for (i = 0; i < PLOVER_N; i++) {
        n += ((__int128_t) z1[i]) * ((__int128_t) z1[i]);
        n += ((__int128_t) sig->z[i]) * ((__int128_t) sig->z[i]);
        n += sig->c1[i] * sig->c1[i];
    }

    return n < PLOVER_BETASQ;
}

//  === plover_core_keygen ===
//  Generate a public-secret keypair ("pk", "sk").

void plover_core_keygen(plover_pk_t *pk, plover_sk_t *sk)
{
    int i, j;
    int64_t a[PLOVER_N];
    int64_t mb[PLOVER_D][PLOVER_N], me[PLOVER_D][PLOVER_N];
    mask_random_t mrg;

    //  intialize the mask random generator
    mask_random_init(&mrg);

    //  --- 1.  seed <- {0,1}^kappa
    randombytes(pk->a_seed, PLOVER_AS_SZ);

    //  --- 2.  a := ExpandA(seed)
    expand_a(a, pk->a_seed);

    //  --- 3.  [[(e, s)]] <- AddRepNoise([[(e, s)]], ut, rep)
    add_rep_noise(me, 0, PLOVER_UT, &mrg);
    add_rep_noise(sk->s, 0, PLOVER_UT, &mrg);
    for (j = 0; j < PLOVER_D; j++) {
        polyr_fntt(sk->s[j]);
    }

    //  --- 4.  [[b]] := A * [[s]] + [[e]]
    for (j = 0; j < PLOVER_D; j++) {
        polyr_ntt_cmul(mb[j], sk->s[j], a);
        polyr_intt(mb[j]);
    }
    for (j = 0; j < PLOVER_D; j++) {
        polyr_addq(mb[j], mb[j], me[j]);
    }

    //  --- 7.  b := Unmask([[b]])
    plover_decode(pk->b, mb, &mrg);

    //  --- 8.  b := round( eta - b )_q->q_t
    for (i = 0; i < PLOVER_N; i++) {
        if (i == 0) {
            pk->b[i] = mont64_sub(pk->b[i], 1L << PLOVER_LOGD);
        }
        pk->b[i] = mont64_cadd(-pk->b[i], PLOVER_Q);
    }

    round_shift_r(pk->b, PLOVER_QT, PLOVER_NUT);
    polyr_nonneg(pk->b, pk->b, PLOVER_Q);

    //  --- 9.  return ( (vk := seed, b'), sk:= (vk, [[s]]) )
    memcpy(&sk->pk, pk, sizeof(plover_pk_t));
}

//  === plover_core_sign ===
//  Create a detached signature "sig" for digest "mu" using secret key "sk".

void plover_core_sign(plover_sig_t *sig, plover_sk_t *sk, const uint8_t *m, 
                    size_t m_sz)
{
    int i, j;
    int64_t a[PLOVER_N];
    int64_t mp[2][PLOVER_D][PLOVER_N];
    int64_t mw[PLOVER_D][PLOVER_N];
    int64_t u[PLOVER_N], c[PLOVER_N], z1[PLOVER_N];
    int64_t w[PLOVER_N], c1_ntt[PLOVER_N];
    int64_t mz1[PLOVER_D][PLOVER_N], mz2[PLOVER_D][PLOVER_N];
    bool rsp = false;
    mask_random_t mrg;

    //  intialize the mask random generator
    mask_random_init(&mrg);

    //  --- 1.  (vk, [[s]]) := [[sk]], (seed, t) := vk      [ caller ]

    //  --- 3.  a := ExpandA(seed)
    expand_a(a, sk->pk.a_seed);

    do {

        //  --- 1.  salt <- {0,1}^{2*kappa}
        randombytes(sig->salt, PLOVER_SALT_SZ);

        xof_msg_hash(u, sig->salt, sk->pk.tr, m, m_sz);

        for (i = 0; i < 2; i++) {
            //  --- 4.  [[p]] <- AddRepNoise([[p]], uw, rep)
            add_rep_noise(mp[i], i, PLOVER_UW, &mrg);
        }

        //  --- 5. [[w]]<- a*[[p_2]] + [[p_1]]
        for (j = 0; j < PLOVER_D; j++) {
            // Transform p[1] into NTT domain
            polyr_fntt(mp[1][j]);

            polyr_ntt_cmul(mw[j], a, mp[1][j]);
            polyr_intt(mw[j]);
            polyr_addq(mw[j], mw[j], mp[0][j]);
        }
        plover_refresh(mw, &mrg);

        // --- 6. w <- Unmask([[w]])
        plover_decode(w, mw, &mrg);

        //  --- 7. c <- u - w
        polyr_subq(c, u, w);

        //  --- 8. (c1, c2) <- Decompose(c)
        int64_t maxc1 = (((PLOVER_Q - 1) >> PLOVER_LOGD) + 1L) >> 1;

        int64_t mid1 = maxc1;
        int64_t mid2 = 1L << (PLOVER_LOGD - 1);
        int64_t mid = (mid1 << PLOVER_LOGD) + mid2;

        for (i = 0; i < PLOVER_N; i++) {
            // First add mid to u to center c_1 and c_2 
            sig->c1[i] = mont64_csub(mont64_add(c[i], mid), PLOVER_Q);

            // compute c_1
            sig->c1[i] = mont64_sub(sig->c1[i] >> PLOVER_LOGD, mid1);
            sig->c1[i] = mont64_cadd(sig->c1[i], PLOVER_Q);
        }

        polyr_copy(c1_ntt, sig->c1);
        polyr_fntt(c1_ntt);

        //  --- 9. [[s]] <- Refresh([[s]])
        plover_ntt_refresh(sk->s, &mrg);

        //  --- 10. [[z_2]] <- c_1*[[s]] + [[p_2]]
        for (i = 0; i < PLOVER_D; i++) {
            // Add the factor R^-1 introduced by the NTT multiplication
            polyr_ntt_smul(mp[1][i], mp[1][i], 1);

            polyr_ntt_mula(mz2[i], c1_ntt, sk->s[i], mp[1][i]);
        }

        //  --- 11. z_2 <- Unmask([[z_2]])
        plover_ntt_decode(sig->z, mz2, &mrg);
        polyr_intt(sig->z);

        polyr_center(sig->z, sig->z, PLOVER_Q);
        polyr_center(sig->c1, sig->c1, PLOVER_Q);

        //  --- 12. z_3 <- c_1
        //  --- 13. z_1' <- u - a*z_2 - b*z_3
        //  --- 14. if (z_1', z_2, z_3) >= B, restart
        rsp = plover_check_bounds(u, a, sk->pk.b, sig);
    } while (!rsp);

    //  --- 16. return sig                                  [caller]
}

//  === plover_core_verify ===
//  Verify that the signature "sig" is valid for digest "mu".
//  Returns true iff signature is valid, false if not valid.
bool plover_core_verify(const plover_sig_t *sig, const plover_pk_t *pk, const uint8_t *m,
                      size_t m_sz)
{
    // int i, j;
    int64_t a[PLOVER_N], u[PLOVER_N];

    //  --- 1.  a := ExpandA(seed)
    expand_a(a, pk->a_seed);

    //  --- 2.  u := H(msg, salt, vk)
    xof_msg_hash(u, sig->salt, pk->tr, m, m_sz);


    //  --- 3. z_1' <- u - a*z_2 - b*z_3
    //  --- 4. if (z_1', z_2, z_3) <= B, accept, else reject
    if (!plover_check_bounds(u, a, pk->b, sig)) {
        return false;
    }

    return true;
}
