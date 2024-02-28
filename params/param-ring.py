"""
Security script for PLOVER-RLWE.
"""
from math import sqrt, log, factorial, pi, exp
import estimator
# For debugging purposes
import sys
import os
import contextlib
import copy
import pprint
from dataclasses import dataclass, field


SIGSK_DEFAULT = (1 << 10)
MIN_QUERIES_DEFAULT = (1 << 35)
MAX_QUERIES_DEFAULT = MIN_QUERIES_DEFAULT
SLACK_B2 = 1.2
MIN_SIGMARED = 0.3


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


def shortest_forgery_fixed_blocksize(n, q, blocksize):
    """
    For a fixed blocksize, compute the shortest forgery an adversary can compute via BKZ + Babai.
    See equation (6) in [CPSWX20].
    """
    delta = delta_func(blocksize)
    L = {}
    for x in range(n):
        y = 3 * n - x
        L[x] = (q ** (n / y)) * (delta ** y)
    m = min(L, key=L.get)
    assert (L[m] == min(L.values()))
    return m, L[m]


def find_minimal_blocksize(n, q, bound):
    """
    Find the minimal BKZ blocksize that allows a successful forgery.
    Remember that the success probability of BKZ is a increasing function of the BKZ blocksize.
    We increase te blocksize the success probability becomes non-negligible.
    """
    # Python has a weird numerical behaviour on lower values
    blocksize = 50
    while(blocksize < 2000):
        blocksize += 1
        m, shortest = shortest_forgery_fixed_blocksize(n, q, blocksize)
        if shortest < bound:
            return m, blocksize
    print("blocksize >= 2000")
    return 0, blocksize


def convert_blocksize_to_security(bkz_blocksize):
    """
    Compute the classical and quantum hardness from the BKZ blocksize
    - Classical: use [BDGL16]
    - Quantum: use [CL21]
    """
    bitsec_classical = round(bkz_blocksize * 0.292, 1)
    bitsec_quantum = round(bkz_blocksize * 0.250, 1)
    return bitsec_classical, bitsec_quantum


@supress_stdout
def LWE_estimate_rough_silent(LWE_instance):
    return estimator.LWE.estimate.rough(LWE_instance)


def compute_RLWE_hardness(n, q, sig_reduc):
    LWEGaussian = estimator.lwe_parameters.NoiseDistribution.DiscreteGaussian
    # Generate and estimate hardness of LWE instance
    LWE_instance = estimator.lwe_parameters.LWEParameters(
        n=n, q=q,
        Xs=LWEGaussian(sig_reduc),
        Xe=LWEGaussian(sig_reduc),
        m=n, tag="",
    )
    LWE_hardness = LWE_estimate_rough_silent(LWE_instance)
    blocksize = min(LWE_hardness[key]["beta"] for key in LWE_hardness)
    # Apply dimensions for free
    blocksize -= dimensionsforfree(blocksize)
    classic, quantum = convert_blocksize_to_security(blocksize)
    return classic


def log2(x):
    """Logarithm in base 2."""
    if (x == 0):
        return -1
    y = log(x, 2)
    if y > sys.float_info.max:
        raise ValueError
    else:
        return y


@dataclass
class CoreParam:
    """
    Dataclass for storing core parameter sets.
    They are used to compte other parameters.
    """
    bit_target: int     # Target bit-security
    n: int              # Dimension of the ring
    q: int              # Modulus
    div: int            # Divider
    sigpert: float      # Standard deviation for pert
    sigsk: float = field(default=MIN_QUERIES_DEFAULT)
    bitdrop: int = field(default=0)
    # Standard deviation for sk (it left empty, computed at initialization)
    min_queries: int = field(default=MIN_QUERIES_DEFAULT)
    # Minimum number of queries (can be left empty, default: MIN_QUERIES_DEFAULT)
    max_queries: int = field(default=MAX_QUERIES_DEFAULT)
    # Maximum number of queries (it left empty, computed at initialization)
    max_queries_above_smooth: int = field(default=MAX_QUERIES_DEFAULT)
    # Maximum number of queries (it left empty, computed at initialization)


@dataclass
class HardnessForge:
    """
    Dataclass for storing metrics related to forgery
    """
    B2: int         # B2 Bound
    BRSIS: int      # RSIS bound
    log2B2: float   # Bound (log)
    forgotten: int  # Forgotten vectors
    blocksize: int  # BKZ blocksize
    classic: int    # Classic security
    quantum: int    # Quantum security


@dataclass
class Sizes:
    """
    Dataclass for storing the size (in bytes) of sig and vk
    """
    vk: int
    sig: int


def printmydataclass(label, obj, log=False):
    """
    Function for pretty-printing dataclass objects
    """
    print("\n" + label)
    xlen = max(len(x) for x in obj.__dataclass_fields__) + 10
    ylen = max(len(str(obj.__getattribute__(x))) for x in obj.__dataclass_fields__) + 10
    if log:
        for x in obj.__dataclass_fields__:
            print("{x} = {y} (log = {z:.2f})".format(
                x=x.ljust(xlen),
                y=str(obj.__getattribute__(x)).ljust(ylen),
                z=log2(obj.__getattribute__(x))
            ))
    else:
        for x in obj.__dataclass_fields__:
            print("{x} = {y}".format(
                x=x.ljust(xlen),
                y=str(obj.__getattribute__(x))
            ))


class Param:
    """
    This class stores all parameters and metrics relative to an instance of PLOVER-RLWE.
    """

    def __init__(self, core, queries=True):
        """
        Initialize a Param object.
        """
        verbose = True
        verbosehardness = True
        self.core = copy.copy(core)
        self.complete_params(verbose=verbose)
        self.compute_forge_hardness(verbosehardness=verbosehardness)
        if (queries is True):
            self.compute_queries(verbose=verbose, verbosehardness=verbosehardness)
        self.compute_sizes()

    def complete_params(self, verbose=True):
        """
        Complete the core parameter set by computing B2 and max_queries
        """
        n = self.core.n
        q = self.core.q
        div = self.core.div
        sigpert = self.core.sigpert

        # Computing and storing bitdrop
        self.core.bitdrop = int(log2((div * sigpert) / (q * sqrt(n)))) - 5

        # Computing and storing sigsk
        self.core.sigsk = 1 << (int(log2((div / q) * sqrt(6 / n) * sigpert)))
        sigsk = self.core.sigsk

        # Computing and storing B2
        # We split B2 in three additive components for easier parameter selection
        B2_pert = int(2 * sigpert * sigpert)
        B2_div = int((div ** 2) / 12)
        B2_sk = int((n / 6) * ((q * sigsk / div) ** 2))
        B2_round = int(n * (2**(2*self.core.bitdrop) / 12) * (q/div)**2 / 12)
        print("B2_round", B2_round)
        if (verbose is True):
            print("B2_pert".ljust(10), B2_pert)
            print("B2_div".ljust(10), B2_div)
            print("B2_sk".ljust(10), B2_sk)
        self.B2 = int( SLACK_B2 * sqrt( n * (B2_pert + B2_div + B2_sk + B2_round) ) )

    def compute_forge_hardness(self, verbose=False, verbosehardness=False):
        """
        Computing the hardness of forgery when cast as an inhomonegeous Ring-SIS instance.
        Based on the state-of-the-art which is https://ia.cr/2019/1456.
        """
        bitdrop = self.core.bitdrop
        B2, log2B2 = self.B2, int(log2(self.B2))
        BRSIS = B2 + int((n ** (3/2)) * q * (2 ** (bitdrop - 2)) / div)

        forgotten, blocksize = find_minimal_blocksize(n, q, B2)
        classic, quantum = convert_blocksize_to_security(blocksize)
        if (classic < self.core.bit_target) and (verbosehardness is True):
            print("Forgery hardness under threshold (Forgery hardness is {classic} bits)".format(classic=classic))
        self.forge = HardnessForge(
            B2=B2,
            BRSIS=BRSIS,
            log2B2=log2B2,
            forgotten=forgotten,
            blocksize=blocksize,
            classic=classic,
            quantum=quantum
        )

    def compute_queries(self, verbose=False, verbosehardness=False):
        """
        Compute the hardness of key-recovery based on the Hint-RLWE -> RLWE reduction (https://ia.cr/2023/623).
        As this hardness is a function of the number of queries, we first set the number to be small,
        then increase it until the bit-security of the RLWE instance is barely higher than bit_target.
        """
        n = self.core.n
        q = self.core.q
        div = self.core.div
        sigsk = self.core.sigsk
        sigpert_PROBED = self.core.sigpert / sqrt(2)
        min_queries = self.core.min_queries
        queries = self.core.min_queries
        bitsec = self.core.bit_target
        eps = 2 ** (- bitsec)
        smoothing_parameter = float(sqrt(log( 4 * n * (1 + 1 / eps)) / 2) / pi)
        if verbose:
            print("Hint-RLWE -> RLWE reduction: (number of queries, bit-security)")
        while(1):
            # Compute sig_reduc
            B_HRLWE = queries * (n / 12) * ((q / div) ** 2)
            sig_reduc = (2 * ( (sigsk ** (- 2)) + B_HRLWE * (sigpert_PROBED ** (- 2)) )) ** (- 1 / 2)
            if (sig_reduc < smoothing_parameter) and (verbosehardness is True):
                print("Hint-RLWE reduction has sigma under smoothing threshold (sigma = {sig_reduc:.3f} < {smooth:.3f} = smoothing)".format(
                    sig_reduc=sig_reduc,
                    smooth=smoothing_parameter)
                )

            # Generate and estimate hardness of LWE instance
            classic = compute_RLWE_hardness(n, q, sig_reduc)

            if verbose:
                print("(2^{x}, {bitsec})".format(x=int(log2(queries)), bitsec=classic))

            # Setting queries
            if (classic >= self.core.bit_target) and (sig_reduc > smoothing_parameter):
                self.core.max_queries_above_smooth = queries
            if (classic >= self.core.bit_target) and (sig_reduc > MIN_SIGMARED):
                self.core.max_queries = queries

            # Break loop
            if (classic < self.core.bit_target) or (sig_reduc < MIN_SIGMARED):
                break

            # Increase the number of queries
            queries <<= 1

        if verbose:
            print("Done!")
        assert(self.core.max_queries > self.core.min_queries)


    def compute_sizes(self):
        """
        Computing the average size of the verification key and the signature.
        """
        bit = self.core.bit_target
        bitdrop = self.core.bitdrop
        n = self.core.n
        q = self.core.q
        div = self.core.div
        sigpert = self.core.sigpert

        # Verification key (vk)
        vk = int(bit / 8 + n * int(log2(q) - bitdrop) / 8)

        # Signature (sig)
        # We estimate the variance of z_2 and z_3
        z2_var = (div ** 2) / 12
        z3_var = sigpert ** 2
        # We estimate the bytesize of z_2 (or z_3) as (n / 8) times the coefficient-wise entropy.
        # For a Gaussian of variance V, its entropy is (1 + log(2 * pi * V)) / 2
        z2_avg = n * (1 + log(2 * pi * z2_var)) / 16
        z3_avg = n * (1 + log(2 * pi * z3_var)) / 16
        sig = int((2 * bit) + z2_avg + z3_avg)

        # Storing the sizes
        self.sizes = Sizes(vk = vk, sig = sig)

    def __repr__(self):
        # Print parameters
        printmydataclass("Core Parameters", self.core, log=True)
        printmydataclass("Hardness of forgery", self.forge)
        printmydataclass("Sizes (indicative)", self.sizes)
        return ""

D = {}
for lq in range(35, 40):
    q = (1 << lq)
    n = 2048
    sigpert = 1 << (lq - 5)
    div = 1 << (lq - 4)
    while(1):
        PLOVER_RLWE_128_core_small = CoreParam(bit_target=128, n=n, q=q, div=div, sigpert=sigpert)
        PLOVER_RLWE_128 = Param(PLOVER_RLWE_128_core_small, queries=False)
        sigpert >>= 1
        div >>= 1
        if (PLOVER_RLWE_128.forge.classic >= PLOVER_RLWE_128.core.bit_target):
            break
    PLOVER_RLWE_128 = Param(PLOVER_RLWE_128_core_small, queries=True)
    print(PLOVER_RLWE_128)
