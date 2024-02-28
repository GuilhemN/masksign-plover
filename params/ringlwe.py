"""
Scripts for Masked Plover's parameters.
References:
- [LatticeE]    https://github.com/malb/lattice-estimator
- [BDGL16]:     https://ia.cr/2015/1128
- [CL21]:       https://ia.cr/2021/570
- [CPSWX20]:    https://ia.cr/2019/1456
- [Duc18]:      https://ia.cr/2017/99
- [Laa16]:      https://pure.tue.nl/ws/files/14673128/20160216_Laarhoven.pdf
- [MW16]:       https://ia.cr/2015/1123
- [NIST]:       https://csrc.nist.gov/CSRC/media/Projects/Post-Quantum-Cryptography
                /documents/call-for-proposals-final-dec-2016.pdf
- [Pre17]:      https://ia.cr/2017/480
"""
from math import sqrt, log, factorial, pi, exp, ceil, floor
from estimator.estimator import lwe_parameters, LWE
import sys
import os
import contextlib


def supress_stdout(func):
    """
    Silence the stdout of a function in Python without
    trashing sys.stdout and restoring each function call.
    See https://stackoverflow.com/a/28321717
    """
    def wrapper(*a, **ka):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                return func(*a, **ka)
    return wrapper


def delta_func(beta):
    """
    Compute delta_beta as per Heuristic 1 in [CPSWX20].
    """
    rep = ((pi * beta) ** (1 / beta)) * beta / (2 * pi * exp(1))
    rep = rep ** (1 / (2 * (beta - 1)))
    return rep


def dimensionsforfree(B):
    """
    Number of "dimensions for free", called d in [Duc18].
    """
    return round(B * log(4 / 3) / log(B / (2 * pi * exp(1))))


def shortest_forgery_fixed_blocksize(k, ell, n, q, blocksize):
    """
    For a fixed blocksize, compute the shortest forgery an adversary
    can compute via BKZ + Babai. See equation (6) in [CPSWX20].
    """
    delta = delta_func(blocksize)
    L = {}
    for m in range(ell * n, (k + ell) * n + 1):
        L[m] = (q ** (k * n / m)) * (delta ** m)
        L[m] = (q ** (k * n / m)) * (delta ** m)
    m = min(L, key=L.get)
    assert (L[m] == min(L.values()))
    return m, L[m]


def find_minimal_blocksize(k, ell, n, q, bound):
    """
    Find the minimal BKZ blocksize that allows a successful forgery.
    The success probability of BKZ increases with the BKZ blocksize.
    We increase te blocksize the success probability becomes non-negligible.
    """
    # Python has a weird numerical behaviour on lower values
    blocksize = 50
    while (blocksize < 2000):
        blocksize += 1
        m, shortest = shortest_forgery_fixed_blocksize(k, ell, n, q, blocksize)
        if shortest < bound:
            return m, blocksize
    print("blocksize >= 2000")
    return 0, blocksize
    # raise ValueError("BKZ blocksize (= {b}) is too big?".format(b=blocksize))


def convert_blocksize_to_security(bkz_blocksize):
    """
    Compute the classical and quantum hardness from the BKZ blocksize
    - Classical: use [BDGL16]
    - Quantum: use [CL21]
    """
    bitsec_classical = int(bkz_blocksize * 0.292)
    bitsec_quantum = int(bkz_blocksize * 0.250)
    return bitsec_classical, bitsec_quantum


def log2(x):
    """Logarithm in base 2."""
    y = log(x, 2)
    if y > sys.float_info.max:
        raise ValueError
    else:
        return y

def smooth(eps, n, normalized=True):
    """
    Compute the smoothing parameter eta_epsilon(Z^n).
    - if normalized is True, take the definition from [Pre17,Falcon]
    - if normalized is False, take the definition from [MR07]
    """
    rep = sqrt(log(2 * n * (1 + 1 / eps)) / pi)
    if normalized is True:
        return rep / sqrt(2 * pi)
    else:
        return rep

class Gaussian:
    def __init__(self, sigma):
        self.label = "Gaussian"
        self.sigma = sigma


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


class Param:
    """
    Class for sets of parameters.
    """
    def __init__(self, label, bit_target, q, logd, sigt, sigp,
                 n, k, ell, bitdrop=0, DFF=True):
        """
        Initiates Raccon parameters with:
        - moduli q, qt, qw
        - ring degree n
        - module dimensions k and ell
        - The noise distributions in t and w
        - The rounding in t and w
        - the flag DFF (for dimensions for free), set by default to True
        """
        assert (isinstance(label, str))
        self.label = label
        self.bit_target = bit_target
        self.q = q
        self.logq = int(log2(self.q))
        self.logd = logd
        self.sigt = sigt
        self.sigp = sigp
        # End of diff
        self.n = n
        self.k = k
        self.ell = ell

        self.bitdrop = bitdrop
        self.DFF = DFF
        self.complete_parameters()

    # This decorator silences stdout
    @supress_stdout
    def compute_keyrec_hardness(self):
        """
        Compute the hardness of distinguishing the public key from uniform.
        General comments:
        - The estimation applies the lattice estimator [LatticeE]
          + "dimensions for free" [Ducas18].
        - We only take into account the LWE noise added by additive errors
          (not by rounding).
        - We approximate the noise distribution by a Gaussian of same width
        """
        # Generating an LWE instance for the lattice estimator
        LWEGaussian = lwe_parameters.NoiseDistribution.DiscreteGaussian
        samples = self.k * self.n
        self.LWE_instance = lwe_parameters.LWEParameters(
            n=self.ell * self.n,
            q=self.q,
            Xs=LWEGaussian(self.sigmlwe),
            Xe=LWEGaussian(self.sigmlwe),
            m=samples,
            tag="Direct key recovery",
        )
        print(self.LWE_instance)

        # We only test the "usvp" and "dual_hybrid" attacks.
        # These are usually the most efficient attacks.
        # For the full test battery, remove the ".rough" suffix.
        self.LWE_hardness = LWE.estimate.rough(self.LWE_instance)
        keyrec_blocksize_usvp = self.LWE_hardness["usvp"]["beta"]
        keyrec_blocksize_dual = self.LWE_hardness["dual_hybrid"]["beta"]
        keyrec_blocksize = min(keyrec_blocksize_usvp, keyrec_blocksize_dual)

        # The dimensions for free (DFF) optimization is not activated
        # by default in the lattice estimator, so we incorporate it manually.
        if (self.DFF is True):
            d = dimensionsforfree(keyrec_blocksize)
            keyrec_blocksize -= d

        # Compute bit security from BKZ blocksize
        self.bit_keyrec_c, self.bit_keyrec_q = convert_blocksize_to_security(keyrec_blocksize)

    def compute_norms(self):
        self.maxc1 = (((self.q-1) >> self.logd) + 1) >> 1

        self.beta = 1.1 * sqrt( self.n * (2 * (self.sigp ** 2) + ((self.maxc1/2) ** 2 + (1 << (self.logd-1)) ** 2) / 4 + (1 << self.bitdrop - 1) * self.maxc1) )
        self.SIS_norm = self.beta + n**(3/2) * (1 << (self.bitdrop-1)) * self.maxc1

    
    def compute_hint_mlwe(self):
        """
        Reduce the security of our scheme to Hint-MLWE, and then MLWE using Theorem 1 of https://eprint.iacr.org/2023/623.pdf.
        """

        sigt = self.sigt
        sigp = self.sigp

        eps = 2**(self.bit_target)/4
        minsigmlwe = sqrt(2)*smooth(eps, self.n)
        self.sigmlwe = 3*minsigmlwe # choose a factor of the minimal sig to increase lwe security, but decrease Q_s

        B = sigp**2 * (1/(2*self.sigmlwe**2)-1/sigt**2)
        bound_approx = 1.001 # tailcut on B such that negligible probability to obtain larger values
        self.mlwe_queries = B / (bound_approx * 4 * self.n * (self.n+3) *(self.maxc1/2)**2)

        # check that bound assumption is correct
        d = sqrt(self.bit_target)

        assert(d * sqrt(self.mlwe_queries/2) * ((2*self.n * self.maxc1)**2) * log(2) <= (bound_approx-1) * self.mlwe_queries * 4 * self.n * (self.n+3) * (self.maxc1/2)**2)

        # Check that sig_p is large enough to hide values statistically in the ROM
        ro_queries = 2**(self.bit_target)
        delta = 2**(-self.bit_target) / ro_queries
        epsilon = (-log(delta, 2) + 8*self.n) / (2*self.n * log(self.q, 2))

        self.minsigp = sqrt(2*pi) * log(4*self.n*(1+1/delta)) / pi * sqrt(n) * self.q**(1/2+epsilon)
        assert(self.sigp > self.minsigp)

        # constraint in new security reduction with \vecz uniform to have enough min entropy in signatures
        # assert(1/sqrt(12) * 2 * self.SIS_norm/sqrt(2*n) > self.minsigp)
        self.min_SIS_norm_red = self.minsigp / (1 / sqrt(12) * 2 / sqrt(2*n))

    def minus_log2_f(self, b, c, x):
        rep = - 2 * sqrt(b * c * x) + b + log(x)
        rep /= log(2)
        return rep

    def compute_forgery_hardness(self, verbose=True):
        """
        Compute the minimal BKZ blocksize for a BKZ + Babai forgery,
        as well as the associated bit security (Classical and Quantum).
        This is based on the forgery analysis in [CPSWX20].
        """
        optimal = self.q ** (self.k / (self.k + self.ell))
        if (self.SIS_norm < optimal):
            if (verbose is True):
                print("This choice of parameters is sub-optimal. Try changing\
                    parameters by e.g. increasing the standard deviation \
                    (or ask Rafael).")
                print("SIS norm =", self.SIS_norm)
                print("optimal  =", optimal)
            self.mm, self.forge_blocksize = 0, 9999
        else:
            self.mm, self.forge_blocksize = find_minimal_blocksize(
                self.k, self.ell, self.n, self.q, self.SIS_norm)
        if (self.DFF is True):
            d = dimensionsforfree(self.forge_blocksize)
            self.forge_blocksize -= d

        # Compute bit security from BKZ blocksize
        self.bit_forge_c, self.bit_forge_q = convert_blocksize_to_security(
            self.forge_blocksize)

    def compute_sizes(self):
        self.seed = int(self.bit_target // 8)
        self.vec = ((self.logq - self.bitdrop) * self.n * self.k) // 8
        self.vk = self.seed + self.vec

        self.sig = 40+ceil(1.1 * log(self.SIS_norm / sqrt(self.n), 2) + 1) * self.n // 8 + self.n * ceil(log(self.maxc1, 2)+1) //8

    def complete_parameters(self):
        """
        Compute all parameters from a set of initial parameters.
        In addition, the bit security, maximum number of queries and
        various sizes are computed.
        """
        self.compute_norms()
        self.compute_hint_mlwe()
        self.compute_keyrec_hardness()
        self.compute_forgery_hardness()
        self.compute_sizes()

    def __repr__(self):
        # Print parameters
        rep = ""
        rep += "ploversign_{bit} = Param(label=\"{label}\", \
                \t\ttarget_bitsec={bit}, q={q}, \
                \n\t\tn={n}, k={k}, ell={ell})\n\n".\
            format(
                label=self.label, bit=self.bit_target, n=self.n, 
                q=self.q, k=self.k, ell=self.ell,
            )
        rep += str(self.label) + "\n\n"

        delimiter = "===============\n"
        rep += "Parameters\n"
        rep += delimiter
        rep += "bit target:         {x}\n".format(x=self.bit_target)
        rep += "q:                  2^{x}\n".format(x=self.logq)
        rep += "d:                  2^{x}\n".format(x=self.logd)
        rep += "sigt:               2^{x}\n".format(x=int(log(self.sigt,2)))
        rep += "sigp:               2^{x} (>= min: 2^{y})\n".format(x=log(self.sigp,2), y=log(self.minsigp, 2))
        rep += "n:                  {x}\n".format(x=self.n)
        rep += "k:                  {x}\n".format(x=self.k)
        rep += "ell:                {x}\n".format(x=self.ell)

        # Print the test results
        rep += "\nLWE check\n"
        rep += delimiter
        rep += "sigmlwe:            {x}\n".format(x=self.sigmlwe)
        rep += "Sec (C):            {x}\n".format(x=self.bit_keyrec_c)
        rep += "Sec (Q):            {x}\n".format(x=self.bit_keyrec_q)
        rep += "Pass?               {x}\n".format(x=(self.bit_keyrec_c > self.bit_target))

        rep += "\nForgery check\n"
        rep += delimiter
        # rep += "sqnorm_average:               2^{x}\n".format(x=log2(self.sqnorm_average))
        # rep += "err norm:                     2^{x}\n".format(x=log2(self.norm_err))
        rep += "beta:                         2^{x}\n".format(x=log2(self.beta))
        rep += "SIS_norm:                     2^{x} (>= 2^{y} for security red)\n".format(x=log2(self.SIS_norm), y=log2(self.min_SIS_norm_red))
        sis_dimension = (self.k + self.ell) * self.n
        rep += "SIS_norm / sqrt(dimension):   2^{x}\n".format(x=log2(self.SIS_norm / sqrt(sis_dimension)))
        rep += "mm:         {x}\n".format(x=self.mm)
        rep += "Blocksize:  {x}\n".format(x=self.forge_blocksize)
        rep += "Sec (C):    {x}\n".format(x=self.bit_forge_c)
        rep += "Sec (Q):    {x}\n".format(x=self.bit_forge_q)
        rep += "Pass?       {x}\n".format(x=(self.bit_forge_c > self.bit_target))

        rep += "\nQueries\n"
        rep += delimiter
        rep += "Queries:              2^{x}\n".format(x=log2(int(self.mlwe_queries)))

        rep += "\nSizes\n"
        rep += delimiter
        rep += "vk:                 {x} bytes\n".format(x=self.vk)
        rep += "sig (INDICATIVE):   {x} bytes\n".format(x=self.sig)

        return rep

n = 2048
q = 549824583172097
logd = 41 # log 2 of the divisor used for signing

parameters_128 = Param(
    label="Parameters 128, SumUniforms",
    bit_target=128, q=q, logd=logd,

    # standard deviation of secret key coefficients
    sigt=SumUniforms(u=6, rep=8).sigma,
    # standard deviation of perturbation coefficients
    sigp=SumUniforms(u=42, rep=8).sigma,

    n=n, k=1, ell=1,
    bitdrop=23,
)
print(parameters_128)
