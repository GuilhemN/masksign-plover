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

static void zero_scalar_encoding(int64_t z[PLOVER_D], mask_random_t *mrg)
{
#if PLOVER_D == 1
    (void)mrg;
    z[0] = 0;
#else
    int i, j, d;
    int64_t r;

    //  d = 2
    for (i = 0; i < PLOVER_D; i += 2)
    {
        z[i] = mask_random_scalar(mrg, i);
        z[i + 1] = mont64_cadd(-z[i], PLOVER_Q);
    }

    //  d = 4, 8, ..
    d = 2;
    while (d < PLOVER_D)
    {
        for (i = 0; i < PLOVER_D; i += 2 * d)
        {
            for (j = i; j < i + d; j++)
            {
                r = mask_random_scalar(mrg, j);
                z[j] = mont64_csub(z[j] + r, PLOVER_Q);
                z[j + d] = mont64_cadd(z[j + d] - r, PLOVER_Q);
            }
        }
        d <<= 1;
    }
#endif
}

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
    (void)x;
    (void)mrg;
#else
    int i;
    int64_t z[PLOVER_D][PLOVER_N];

    //  --- 1.  [[z]] <- ZeroEncoding(d)
    zero_encoding(z, mrg);

    //  --- 2.  return [[x]]' := [[x]] + [[z]]
    for (i = 0; i < PLOVER_D; i++)
    {
        polyr_addq(x[i], x[i], z[i]);
    }
#endif
}

static void plover_scalar_refresh(int64_t x[PLOVER_D], mask_random_t *mrg)
{
#if PLOVER_D == 1
    (void)x;
    (void)mrg;
#else
    int i;
    int64_t z[PLOVER_D];

    //  --- 1.  [[z]] <- ZeroEncoding(d)
    zero_scalar_encoding(z, mrg);

    //  --- 2.  return [[x]]' := [[x]] + [[z]]
    for (i = 0; i < PLOVER_D; i++)
    {
        x[i] = mont64_csub(x[i] + z[i], PLOVER_Q);
    }
#endif
}

//  Decode(): Collapse shares

static void plover_decode(int64_t r[PLOVER_N], int64_t m[PLOVER_D][PLOVER_N], mask_random_t *mrg)
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

//  Decode(): Collapse shares (possibly split CRT arithmetic)

static void plover_ntt_decode(int64_t r[PLOVER_N], int64_t m[PLOVER_D][PLOVER_N], mask_random_t *mrg)
{
    plover_ntt_refresh(m, mrg);

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

static void plover_scalar_masked_mul(int64_t p[PLOVER_D], const int64_t x[PLOVER_D], const int64_t y[PLOVER_D], mask_random_t *mrg)
{
    int i, j;
    int64_t s[PLOVER_D], ret[PLOVER_D];
    for (i = 0; i < PLOVER_D; i++)
    {
        ret[i] = mont64_cadd(mont64_mulq(x[i], y[i]), PLOVER_Q);
    }
    for (i = 0; i < PLOVER_D; i++)
    {
        zero_scalar_encoding(s, mrg); // Draw a random masking polynomial
        for (j = i + 1; j < PLOVER_D; j++)
        {
            ret[j] = mont64_cadd(ret[j] - s[j], PLOVER_Q);

            s[j] = mont64_csub(mont64_cadd(mont64_mulq(x[i], y[j]), PLOVER_Q) + s[j],
                               PLOVER_Q);
            s[j] = mont64_csub(mont64_cadd(mont64_mulq(x[j], y[i]), PLOVER_Q) + s[j],
                               PLOVER_Q);

            ret[i] = mont64_csub(ret[i] + s[j], PLOVER_Q);
        }
    }

    for (i = 0; i < PLOVER_D; i++)
    {
        p[i] = ret[i];
    }
}

//  exponentiation of a masked scalar
//  Non constant-time with regard to e: assumes that e is public

// Returns a result in Montgomery space

static void plover_scalar_masked_pow(int64_t y[PLOVER_D], const int64_t x[PLOVER_D], int64_t e, mask_random_t *mrg)
{
    int j;
    int64_t x2[PLOVER_D], tmp[PLOVER_D], tmpx2[PLOVER_D];
    // x in Montgomery domain
    for (j = 0; j < PLOVER_D; j++)
    {
        x2[j] = mont64_cadd(mont64_mulq(x[j], MONT_RR), PLOVER_Q);
    }

    // Set y to the polynomial with all coefficients equal to 1 in Montgomery domain
    y[0] = MONT_R;
    for (j = 1; j < PLOVER_D; j++)
    {
        y[j] = 0;
    }

    while (e > 0)
    {
        if ((e & 1) == 1)
        {
            plover_scalar_masked_mul(tmp, y, x2, mrg); // [[y]] <- [[y]] * [[x]]
            for (j = 0; j < PLOVER_D; j++) {
                y[j] = tmp[j];
            }
        }

        for (j = 0; j < PLOVER_D; j++) {
            tmpx2[j] = x2[j];
        }
        plover_scalar_refresh(tmpx2, mrg);
        plover_scalar_masked_mul(tmp, x2, tmpx2, mrg); // [[x]] <- [[x]] * [[x]]
        for (j = 0; j < PLOVER_D; j++) {
            x2[j] = tmp[j];
        }
        e >>= 1;
    }
}

//  coordinate wise multiplication of two masked vectors
//  p must be different from x and y

static void plover_masked_mul(int64_t p[PLOVER_D][PLOVER_N], int64_t x[PLOVER_D][PLOVER_N], int64_t y[PLOVER_D][PLOVER_N], mask_random_t *mrg)
{
    XASSERT(p != x && p != y); // not supported by the code

    int i, j;
    int64_t s[PLOVER_N];
    for (i = 0; i < PLOVER_D; i++) {
        polyr_ntt_cmul(p[i], x[i], y[i]);
    }
    for (i = 0; i < PLOVER_D; i++) {
        for (j = i + 1; j < PLOVER_D; j++) {
            mask_random_poly(mrg, s, i); // Draw a random masking polynomial

            polyr_subq(p[j], p[j], s);

            polyr_ntt_mula(s, x[i], y[j], s);
            polyr_ntt_mula(s, x[j], y[i], s);

            polyr_addq(p[i], p[i], s);
        }
    }
}


//  coordinate wise exponentation of a vector
//  Non constant-time with regard to e: assumes that e is public

// Returns a result in Montgomery space

static void plover_masked_pow(int64_t y[PLOVER_D][PLOVER_N], const int64_t x[PLOVER_D][PLOVER_N], int64_t e, mask_random_t *mrg)
{
    int i, j;
    int64_t tmp[PLOVER_D][PLOVER_N], x2[PLOVER_D][PLOVER_N], tmpx2[PLOVER_D][PLOVER_N];
    // x in Montgomery domain
    for (j = 0; j < PLOVER_D; j++) {
        polyr_ntt_smul(x2[j], x[j], MONT_RR);
    }

    // Set y to the polynomial with all coefficients equal to 1 in Montgomery domain
    for (i = 0; i < PLOVER_N; i++) {
        y[0][i] = MONT_R;
        for (j = 1; j < PLOVER_D; j++) {
            y[j][i] = 0;
        }
    }

    while (e > 0) {
        if ((e & 1) == 1) {
            plover_masked_mul(tmp, y, x2, mrg); // [[y]] <- [[y]] * [[x]]
            for (j = 0; j < PLOVER_D; j++) {
                polyr_copy(y[j], tmp[j]);
            }
        }

        for (j = 0; j < PLOVER_D; j++) {
            polyr_copy(tmpx2[j], x2[j]);
        }
        plover_refresh(tmpx2, mrg);
        plover_masked_mul(tmp, x2, tmpx2, mrg); // [[x]] <- [[x]] * [[x]]
        for (j = 0; j < PLOVER_D; j++) {
            polyr_copy(x2[j], tmp[j]);
        }

        e >>= 1;
    }
}

static int64_t pseudo_inverse(int64_t mg_inv[PLOVER_D][PLOVER_N], int64_t mg_ntt[PLOVER_D][PLOVER_N], mask_random_t *mrg)
{
    //  --- 5. [[g_\times]] <- PseudoInverse([[g]])
    // Compute g_inv <- g^{-1} = g^{phi(q)-1}

    int i, j;
    int64_t me_check[PLOVER_D], e_check;
    int64_t prod_g[PLOVER_N][PLOVER_D]; // Compute partial products of coordinates of g
    // prod_g[i] = prod_{j <= i} g[j] * R^(-i)
    // First coordinate
    for (i = 0; i < PLOVER_D; i++)
    {
        prod_g[0][i] = mg_ntt[i][0];
    }
    // other partial products
    for (i = 1; i < PLOVER_N; i++)
    {
        int64_t gi[PLOVER_D];
        for (j = 0; j < PLOVER_D; j++)
        {
            gi[j] = mg_ntt[j][i];
        }
        plover_scalar_masked_mul(prod_g[i], prod_g[i - 1], gi, mrg);
    }

    // Compute inverse of product of all coordinates
    // prod_g_inv = (prod g[j]) * R^n
    int64_t prod_g_inv[PLOVER_D];
    plover_scalar_masked_pow(prod_g_inv, prod_g[PLOVER_N - 1], PLOVER_PHIQ - 1, mrg);

    // Recover inverse of each coordinate
    int64_t tmp_mg_inv[PLOVER_N][PLOVER_D];
    int64_t prod_g_end[PLOVER_D]; // Partial product of g coordinates, starting from the end
    prod_g_end[0] = MONT_R;
    for (j = 1; j < PLOVER_D; j++)
    {
        prod_g_end[j] = 0;
    }

    plover_refresh(mg_ntt, mrg); // Refresh before re-use

    for (i = PLOVER_N - 1; i >= 0; i--)
    {
        // Recover inverse of coordinate i
        // mginv2[i] <- (prod g[j])^(-1) * R^{n-1}
        memcpy(tmp_mg_inv[i], prod_g_inv, sizeof(prod_g_inv));

        if (i > 0) {
            // mginv2[i] = (prod_{j >= i})^(-1) * R^(n-1-i)
            plover_scalar_masked_mul(tmp_mg_inv[i], tmp_mg_inv[i], prod_g[i - 1], mrg);
        }

        plover_scalar_masked_mul(tmp_mg_inv[i], tmp_mg_inv[i], prod_g_end, mrg);

        int64_t gi[PLOVER_D];
        for (j = 0; j < PLOVER_D; j++)
        {
            gi[j] = mg_ntt[j][i];
        }
        plover_scalar_refresh(prod_g_end, mrg);
        plover_scalar_masked_mul(prod_g_end, prod_g_end, gi, mrg);
    }


    for (i = 0; i < PLOVER_N; i++) {
        for (j = 0; j < PLOVER_D; j++) {
            mg_inv[j][i] = tmp_mg_inv[i][j];
        }
    }

    plover_scalar_refresh(prod_g[PLOVER_N - 1], mrg); // Refresh before re-use
    plover_scalar_masked_mul(me_check, prod_g_inv, prod_g[PLOVER_N-1], mrg);
    
    // Unmasked e, if g is inversible e=1, otherwise e=0
    e_check = 0;
    for (j = 0; j < PLOVER_D; j++) {
        e_check = mont64_csub(e_check + me_check[j], PLOVER_Q);
    }

    mask_random_t mrg2;
    mask_random_init(&mrg2);

    return e_check;
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

//  CheckBounds(sig) -> {OK or FAIL}

static bool plover_check_bounds( const int64_t u[PLOVER_N], const int64_t h_ntt[PLOVER_N], const int64_t z2[PLOVER_N])
{
    int i;
    int64_t z1[PLOVER_N], z2_ntt[PLOVER_N];

    // Recompute z_1
    polyr_copy(z2_ntt, z2);
    polyr_fntt(z2_ntt);

    // Compute z1 <- u - h * z_2 = u' + Ap - h z_2 = u' + p_1 - h c_1 * g
    //                      = u' + p_1 - c_1 (eta - f)
    polyr_ntt_cmul(z1, h_ntt, z2_ntt);
    polyr_intt(z1);

    polyr_subq(z1, u, z1);

    // Check norm on vector z
    polyr_center(z1, z1, PLOVER_Q);
    __int128_t n = 0;
    for (i = 0; i < PLOVER_N; i++) {
        n += ((__int128_t) z1[i]) * ((__int128_t) z1[i]);
        n += ((__int128_t) z2[i]) * ((__int128_t) z2[i]);
    }

    return n < PLOVER_BETASQ;
}


//  === plover_core_keygen ===
//  Generate a public-secret keypair ("pk", "sk").

void plover_core_keygen(plover_pk_t *pk, plover_sk_t *sk)
{
    int i, j;
    int64_t mh[PLOVER_D][PLOVER_N], mf[PLOVER_D][PLOVER_N], mg_inv[PLOVER_D][PLOVER_N];
    int64_t me_check[PLOVER_D], e_check;
    mask_random_t mrg;

    //  intialize the mask random generator
    mask_random_init(&mrg);

    //  --- 2.  [[f]] <- AddRepNoise([[f]], ut, rep)
    add_rep_noise(mf, 0, PLOVER_UT, &mrg); // compute f' <- f - eta
    mf[0][0] = mont64_cadd(mont64_sub(mf[0][0], 1L << PLOVER_LOGD), PLOVER_Q);

    for (j = 0; j < PLOVER_D; j++) {
        polyr_fntt(mf[j]);
    }

    //  --- 3.  while(e = 0)
    do {
        //  --- 4.   [[g]] <- AddRepNoise([[g]], ut, rep)
        add_rep_noise(sk->g_ntt, 0, PLOVER_UT, &mrg);
        for (j = 0; j < PLOVER_D; j++) {
            polyr_fntt(sk->g_ntt[j]);
        }

        e_check = pseudo_inverse(mg_inv, sk->g_ntt, &mrg);
    } while (e_check == 0);

    //  --- 6.   [[h]] <- (divider-f)/g
    plover_masked_mul(mh, mf, mg_inv, &mrg);
    plover_ntt_decode(pk->h, mh, &mrg);

    // h <- -h
    polyr_ntt_smul(pk->h, pk->h, 1);
    polyr_intt(pk->h);
    polyr_negm(pk->h, pk->h, PLOVER_Q);

    //  --- 9.  return ( (vk := seed, h), sk:= (vk, [[g]]) )
    memcpy(&sk->pk, pk, sizeof(plover_pk_t));
}

//  === plover_core_sign ===
//  Create a detached signature "sig" for digest "mu" using secret key "sk".

void plover_core_sign(plover_sig_t *sig, plover_sk_t *sk, const uint8_t *m,
                    size_t m_sz)
{
    int i, j;
    int64_t h_ntt[PLOVER_N];
    int64_t mp[2][PLOVER_D][PLOVER_N];
    int64_t mw[PLOVER_D][PLOVER_N];
    int64_t w[PLOVER_N];
    int64_t u[PLOVER_N], c[PLOVER_N];
    int64_t mz2[PLOVER_D][PLOVER_N];
    int64_t c1[PLOVER_N];
    bool rsp = false;
    mask_random_t mrg;

    //  intialize the mask random generator
    mask_random_init(&mrg);

    //  --- 1.  (vk, [[s]]) := [[sk]], (seed, t) := vk      [ caller ]

    polyr_copy(h_ntt, sk->pk.h);
    polyr_fntt(h_ntt);

    do {

        //  --- 1.  salt <- {0,1}^{2*kappa}
        randombytes(sig->salt, PLOVER_SALT_SZ);

        xof_msg_hash(u, sig->salt, sk->pk.tr, m, m_sz);

        for (i = 0; i < 2; i++) {
            //  --- 3.  [[p]] <- AddRepNoise([[p]], uw, rep)
            add_rep_noise(mp[i], i, PLOVER_UW, &mrg);
        }

        //  --- 4.  [[w]]<- [ 1  h ]*[[p]]
        for (j = 0; j < PLOVER_D; j++) {
            // Transform p[1] into NTT domain
            polyr_fntt(mp[1][j]);

            polyr_ntt_cmul(mw[j], h_ntt, mp[1][j]);
            polyr_intt(mw[j]);
            polyr_addq(mw[j], mw[j], mp[0][j]);
        }

        //  --- 5.  w <- Unmask([[w]])
        plover_decode(w, mw, &mrg);

        //  --- 6. c <- u - w
        polyr_subq(c, u, w);

        //  --- 7. (c_1, c_2) <- Decompose(c)
        int64_t maxc1 = (((PLOVER_Q - 1) >> PLOVER_LOGD) + 1L) >> 1;

        int64_t mid1 = maxc1;
        int64_t mid2 = 1L << (PLOVER_LOGD - 1);
        int64_t mid = (mid1 << PLOVER_LOGD) + mid2;

        for (i = 0; i < PLOVER_N; i++) {
            // First add mid to u to center c_1 and c_2 
            c1[i] = mont64_csub(mont64_add(c[i], mid), PLOVER_Q);

            // compute c_1
            c1[i] = mont64_sub(c1[i] >> PLOVER_LOGD, mid1);
            c1[i] = mont64_cadd(c1[i], PLOVER_Q);
        }

        polyr_fntt(c1);

        //  --- 8.  [[g]] <- Refresh([[g]])
        plover_ntt_refresh(sk->g_ntt, &mrg);

        //  --- 9.  [[z_2]] <- c_1*[[g]] + [[p_2]]
        for (i = 0; i < PLOVER_D; i++) {
            // Add the factor R^-1 introduced by the NTT multiplication
            polyr_ntt_smul(mp[1][i], mp[1][i], 1);

            polyr_ntt_mula(mz2[i], c1, sk->g_ntt[i], mp[1][i]);
        }

        //  --- 9.  z_2 := Unmask([[z_2]])
        plover_ntt_decode(sig->z2, mz2, &mrg);

        polyr_intt(sig->z2);
        polyr_center(sig->z2, sig->z2, PLOVER_Q);

        //  --- 10.  z_1' <- u-h*z_2
        //  --- 11.  if ||(z_1', z_2)|| > B, restart
        rsp = plover_check_bounds(u, h_ntt, sig->z2);
    } while (!rsp);

    //  --- 21. return sig                                  [caller]
}

//  === plover_core_verify ===
//  Verify that the signature "sig" is valid for digest "mu".
//  Returns true iff signature is valid, false if not valid.
bool plover_core_verify(const plover_sig_t *sig, const plover_pk_t *pk, const uint8_t *m,
                      size_t m_sz)
{
    int64_t u[PLOVER_N], h_ntt[PLOVER_N];

    polyr_copy(h_ntt, pk->h);
    polyr_fntt(h_ntt);

    xof_msg_hash(u, sig->salt, pk->tr, m, m_sz);

    //  --- 2.  z_1' <- u-h*z_2
    //  --- 3.  if ||(z_1', z_2)|| > B, restart
    if (!plover_check_bounds(u, h_ntt, sig->z2)) {
        return false;
    }

    return true;
}
