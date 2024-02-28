from sage.all import *

MAX_LOG_N = 13

def valq(candidate):
    (logqi, l, k) = candidate
    return 2**logqi - l*(2**k) + 1

def reprq(candidate):
    (logqi, l, k) = candidate
    return f"2**{logqi} - {l}*2**{k}"

"""
Generate candidates q_i, such that q_i-1 is dividable by 2**k.
"""
def gen_prime_candidates(logqi, k):
    l = 1
    while l <= 2000:
        candidate = (logqi, l, k)
        if valq(candidate) <= 2**(logqi-1):
            break
        yield candidate

        l += 2

"""
Find q a multiple of two primes q_1 and q_2 such that:
    - q_1 and q_2 are dividible by a large power of two
    - log(q_1*q_2) is close to a target
"""
def find_q(logq):
    logq1 = logq // 2
    logq2 = logq - logq1

    cont = True
    for logq1 in range(logq//2, MAX_LOG_N+2-1, -1):
        logq2 = logq - logq1

        # for k in range(logq1-1, MAX_LOG_N+1-1, -1):
        for k in range(MAX_LOG_N+1-1, logq1-1):
            q1, q2 = None, None
            for candidate in gen_prime_candidates(logq1, k):
                if is_prime(valq(candidate)):
                    q1 = candidate
                    break
            for candidate in gen_prime_candidates(logq2, k):
                if is_prime(valq(candidate)) and q1 != candidate:
                    q2 = candidate
                    break
            
            if q1 != None and q2 != None:
                return (q1, q2)

"""
Find h the generator of a group of order 2*n in Z_q.
"""
def find_h(q, n):
    (q1, q2) = q
    q1 = valq(q1)
    q2 = valq(q2)
    q = q1*q2

    F1 = GF(q1)
    F2 = GF(q2)

    # We search the smallest x > 1 generating Z_{q_1}* and Z_{q_2}*
    x = 2
    while True:
        m1 = F1(x).multiplicative_order()
        m2 = F2(x).multiplicative_order()

        if m1 == q1-1 and m2 == q2-1:
            break

        x += 1

    # order of x in Z_q
    m = lcm(m1, m2)

    # We derive an element of 2*n in Z_q
    assert(m % 2*n == 0) # 2*n divides m
    h = pow(x, m//(2*n), q)

    assert((h**n) % q - q == -1)
    assert((h**(2*n))%q == 1)

    return h


qs = {logq: find_q(logq) for logq in range(35, 50)}
assert(all([qs[logq] is not None for logq in qs]))

class SumUniforms:
    """
    Distribution of the sum of rep uniform samples, each having u bits
    """
    def __init__(self, u, rep):
        assert (rep >= 2), "Don't set rep < 2"
        self.label = "SumUniforms"
        self.u = u
        self.rep = rep
        self.sigma = (1 << self.u) * sqrt(rep / 12)

# Fix rep depending on masking order d
reps = {
    1: 8,
    2: 4,
    4: 2,
    8: 4,
    16: 2,
    32: 4,
}

"""
Choose number of bits sampled for the uniforms depending on the masking 
order d, and the target standard deviation sig.
"""
def choose_u(d, logsig, rep):
    # the variance of the sum of uniforms is rep * d * 2^{2u} / 12 = 2^{2logsig}
    # thus u = 1/2 * log(12 * 2^{2logsig} / (rep * d))
    u = 1/2 * log(12 * 2**(2*logsig) / (rep * d)) / log(2)
    return floor(u)

class ConcreteParameters:
    def __init__(self, bitsec, logn, logq, logdivider, logsigpert, logsigsk, nu) -> None:
        self.bitsec = bitsec
        self.logn = logn
        self.n = 2**self.logn
        self.logq = logq
        self.logdivider = logdivider
        self.logsigpert = logsigpert
        self.logsigsk = logsigsk
        self.nu = nu

        self.compute_params()
        self.compute_sizes()

    def compute_params(self):
        assert(self.logq in qs)
        self.q1, self.q2 = qs[self.logq]
        self.qval = valq(self.q1) * valq(self.q2)

        n = 2**self.logn
        self.h = find_h((self.q1, self.q2), n)

        # Compute parameters for AddRepNoise
        self.masking_params = {}
        for d in reps:
            rep = reps[d]
            usk = choose_u(d, logsig=self.logsigsk, rep=rep)
            upert = choose_u(d, logsig=self.logsigpert, rep=rep)                            
            self.masking_params[d] = {
                "rep": rep,
                "usk": usk,
                "upert": upert,
            }

        # compute parameters for NTT
        self.compute_w()

        self.maxc1 = (((self.qval-1) >> self.logdivider) + 1) >> 1

        self.sq_beta = self.n * (2 * (2**(2*self.logsigpert)) + 2**(2*self.logdivider) / 12 + self.qval**2 * self.n / (6 * 2**(2*self.logdivider)) * 2**(2*self.logsigsk))

        print("sq_beta", float(log(self.sq_beta)/log(2)))

        # additional error introduced by the rounding of t
        self.sq_beta += self.n**2 * ((1<<self.nu)*self.maxc1)**2

        self.sq_beta = ceil(self.sq_beta * 1.2**2)

    def compute_w(self):

        def _modexp(x, e, n):
            """(TESTING) Modular exponentiation: Compute x**e (mod n)."""
            y = 1
            while e > 0:
                if e & 1 == 1:
                    y = (y * x) % n
                x = (x * x) % n
                e >>= 1
            return y

        def _bitrev(x, l):
            """(TESTING) Return x with bits 0,1,..(l-1) in reverse order."""
            y = 0
            for i in range(l):
                y |= ((x >> i) & 1) << (l - i - 1)
            return y

        #   g=15 is the smallest generator of both prime fields of composite q.
        #   Reduce to subgroup of order 2*n to obtain "h" in gp-pari;
        #   g   = Mod(15, q)
        #   h   = g^(znorder(g)/(2*n))

        """(TESTING) Re-generate the NTT "tweak" table."""
        q   = self.qval
        lgn = self.logn                 #   log2(n)
        n   = 2**lgn                    #   length of the transform
        h   = self.h                    #   Generates a subgroup of order 2*n
                                        #   obtained with test_params.py
        self.w   = []
        for i in range(n):
            j = _bitrev(i, lgn)
            x = (_modexp(h, j, q)) % q
            self.w.append(x)

    def compute_sizes(self):
        #   nist serialization sizes
        self.pk_sz  =   (self.bitsec // 8 +
                            self.n * (self.logq - self.nu) // 8)
        for d in self.masking_params:
            self.masking_params[d]["sk_sz"] = (self.pk_sz + (d - 1) * (self.bitsec // 8) +
                    (self.n * self.logq) // 8)

        self.sig_sz = 14000
        self.sig_rate = 37

        self.enc_c1_bits = ceil(log(self.maxc1)/log(2))

    def __repr__(self) -> str:
        s = ""
        s += f"Concrete parameters: \n"
        s += f"    n = {2**self.logn}\n"
        s += f"    q = {self.qval} = {valq(self.q1)} * {valq(self.q2)} (log = {float(log(self.qval)/log(2))})\n"
        s += f"    sig_rate = {self.sig_rate}\n"
        s += f"    enc_c1_bits = {self.enc_c1_bits}\n"
        s += f"    sig_sz = {self.sig_sz}\n"


        s += f"\n"
        s += f"Masking parameters: \n"
        for d in self.masking_params:
            p = self.masking_params[d]
            s += f"    d = {d}: rep = {p['rep']}, usk = {p['usk']}, upert = {p['upert']}\n"
        return s
    
    def _repr_w(self):
        w = self.w

        s = ""
        for i in range(0, 2**self.logn, 4):
            s += f"\t{w[i]:15}, {w[i+1]:15}, {w[i+2]:15}, {w[i+3]:15},\n"
        return s

    def _repr_w_for_c(self):
        w = [(c*2**64)%self.qval for c in self.w]

        s = ""
        s += f"\t{w[1]:15}, {w[2]:15}, {w[3]:15},\n"
        for i in range(4, 2**self.logn, 4):
            s += f"\t{w[i]:15}, {w[i+1]:15}, {w[i+2]:15}, {w[i+3]:15},\n"
        return s

    def gen_python(self):
        with open("generated_ring/polyr_params.py", "w") as f:
            f.write(f"# file generated with params/concrete-param-ring.py\n\n")
            f.write(f"PLOVERSIGN_LOGN = {self.logn}\n")
            f.write(f"PLOVERSIGN_N = 2**PLOVERSIGN_LOGN\n")
            f.write(f"PLOVERSIGN_Q1 = {valq(self.q1)}\n")
            f.write(f"PLOVERSIGN_Q2 = {valq(self.q2)}\n")
            f.write(f"PLOVERSIGN_Q = PLOVERSIGN_Q1 * PLOVERSIGN_Q2\n")
            f.write(f"PLOVERSIGN_NI = pow(PLOVERSIGN_N, -1, PLOVERSIGN_Q)   #   n^-1  (mod q)\n")

            f.write(f"\n")
            f.write(f"PLOVERSIGN_H = {self.h}\n")
            f.write(f"PLOVERSIGN_W = [\n{self._repr_w()}]\n")

        with open("generated_ring/ploversign_sets.py", "w") as f:
            f.write(f"# generated with params/gen_concrete.py\n\n")
            for d in self.masking_params:
                p = self.masking_params[d]
                f.write(f"plover_128_{d}  = NIST_PloverSign(  bitsec=128, q=PLOVERSIGN_Q, logdivide={self.logdivider}, nut={self.nu}, rep={p['rep']}, ut={p['usk']}, \n\t\t\t\tup={p['upert']}, n=PLOVERSIGN_N, d={d}, enc_c1_bits={self.enc_c1_bits}, sig_rate={self.sig_rate}, sig_sz={self.sig_sz})\n\n")
            f.write(f"# end generated\n")

    def gen_c(self):
        n = 2**self.logn

        r = (2**64) % self.qval
        rr = (r*r) % self.qval
        ni = (rr * pow(n, -1, self.qval)) % self.qval
        qi = pow(-self.qval, -1, 2**64)

        # approximate sq_beta for the C code
        sq_beta_pow = 0
        sq_beta_high = self.sq_beta
        while sq_beta_high > 2**10:
            sq_beta_high >>= 1
            sq_beta_pow += 1

        with open("generated_ring/mont64.h", "w") as f:
            f.write(f"""// file generated with params/concrete-param-ring.py

#if (PLOVER_N != {n} || PLOVER_Q != {self.qval}l)
#error "Unrecognized polynomial parameters N, Q"
#endif

/*
    n   = {n}
    q1  = {reprq(self.q1)}
    q2  = {reprq(self.q2)}
    q   = q1*q2
    r   = 2^64 % q
    rr  = r^2 % q
    ni  = lift(rr * Mod(n,q)^-1)
    qi  = lift(Mod(-q,2^64)^-1)
*/

//  Montgomery constants. These depend on Q and N
#define MONT_R {r}L
#define MONT_RR {rr}L
#define MONT_NI {ni}L
#define MONT_QI {qi}L        

// end generated      
""")
            
        with open("generated_ring/ntt64.c", "w") as f:
            f.write(f"""// file generated with params/concrete-param-ring.py

static const int64_t plover_w_64[{n-1}] = {{
{self._repr_w_for_c()}}};

// end generated      
""")
            
        with open("generated_ring/param_list.h", "w") as f:
            f.write(f"""// file generated with params/concrete-param-ring.py

#define PLOVER_KAPPA  {self.bitsec}
#define PLOVER_Q      {self.qval}L
#define PLOVER_N      {n}
#define PLOVER_LOGD   {self.logdivider}
#define PLOVER_BETASQ ((__int128_t) {sq_beta_high}) << {sq_beta_pow}

#define PLOVER_SIG_BITS_C1 {self.enc_c1_bits}
#define PLOVER_SIG_RATE    {self.sig_rate}

#define PLOVER_PK_SZ  5136
#define PLOVER_SIG_SZ 14000
""")
            for d in self.masking_params:
                p = self.masking_params[d]
                f.write(f"""
#if defined(PLOVER_128_{d})
#define PLOVER_(s)    PLOVER_128_{d}_##s
#define PLOVER_NAME   "Plover-128-{d}"
#define PLOVER_NUT    {self.nu}
#define PLOVER_REP    {p['rep']}
#define PLOVER_UT     {p['usk']}
#define PLOVER_UW     {p['upert']}
#define PLOVER_D      {d}
#define PLOVER_SK_SZ  {p['sk_sz']}
#endif
""")
            
            f.write("\n// end generated")


p = ConcreteParameters(
    bitsec=128, 
    logn=11, # n = 2048 
    logq=41, 
    logdivider=37, 
    logsigpert=36, 
    logsigsk=27, 
    nu=21
)
print(p)

p.gen_python()
p.gen_c()