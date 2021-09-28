# For licensing see accompanying LICENSE file.
# Copyright (C) 2020 Apple Inc. All Rights Reserved.

import scipy.stats as stats
import math
import numpy as np

# This document contains 4 computations: 2 empirical and 2 theoretical.
# 1. Empirical analysis
# 2. Theoretical analysis

# ========= SUPPORT FUNCTIONS ==========

# This function uses binary search to approximate the smallest eps such that deltacomp will output something smaller than delta (i.e. an algorithm is (eps, delta)-DP)
def binarysearch(deltacomp, delta, num_iterations, epsupper):
    '''
    binary search to find min epsilon such that deltacomp(epsilon)<delta
    deltacomp = function that takes epsilon as input and outputs delta
    num_iterations = number of iterations, accuracy is 2^(-num_iterations)*epsupper
    epsupper = upper bound for epsilon. You should be sure that deltacomp(epsupper)<delta.
    '''
    llim = 0
    rlim = epsupper
    for t in range(num_iterations):
        mideps = (rlim + llim) / 2
        delta_for_mideps = deltacomp(mideps, delta)
        if delta_for_mideps < delta:
            llim = llim
            rlim = mideps
        else:
            llim = mideps
            rlim = rlim
    return rlim

# ================/EXACT EMPIRICAL ANALYSIS WITH STEPS - SAMPLING EMPIRICAL/==============

#This a subroutine in the main algorithm.
def onestep(c, eps, eps0, pminusq):
    '''
    onestep computes the e^(eps)-divergence between p=alpha*Bin(c,0.5)+(1-alpha)*(Bin(c,1/2)+1) and q=alpha*(Bin(c,0.5)+1)+(1-alpha)*Bin(c,1/2), where alpha=e^(eps)/(1+e^(eps))
    if pminusq=True then computes D_(e^eps)(p|q), else computes D_(e^eps)(q|p)
    '''
    alpha = math.exp(eps0) / (math.exp(eps0) + 1)
    effeps = math.log(((math.exp(eps) + 1) * alpha - 1) / ((1 + math.exp(eps)) * alpha - math.exp(eps)))
    if pminusq == True:
        beta = 1 / (math.exp(effeps) + 1)
    else:
        beta = 1 / (math.exp(-effeps) + 1)
    cutoff = beta * (c + 1)
    pconditionedonc = (alpha * stats.binom.cdf(cutoff, c, 0.5) + (1 - alpha) * stats.binom.cdf(cutoff - 1, c, 0.5))
    qconditionedonc = ((1 - alpha) * stats.binom.cdf(cutoff, c, 0.5) + alpha * stats.binom.cdf(cutoff - 1, c, 0.5))
    if pminusq == True:
        return (pconditionedonc - math.exp(eps) * qconditionedonc)
    else:
        return ((1 - qconditionedonc) - math.exp(eps) * (1 - pconditionedonc))


def deltacomp(n, eps0, eps, deltaupper, step, upperbound = True):
    '''
    Let C=Bin(n-1, e^(-eps0)) and A=Bin(c,1/2) and B=Bin(c,1/2)+1 and alpha=e^(eps0)/(e^(eps0)+1)
    p samples from A w.p. alpha and B otherwise
    q samples from B w.p. alpha and A otherwise
    deltacomp attempts to find the smallest delta such P and Q are (eps,delta)-indistinguishable, or outputs deltaupper if P and Q are not (eps, deltaupper)-indistinguishable.
    If upperbound=True then this produces an upper bound on the true delta (except if it exceeds deltaupper), and if upperbound=False then it produces a lower bound.
    '''
    deltap = 0  # this keeps track of int max{0, p(x)-q(x)} dx
    deltaq = 0  # this keeps track of int max{0, q(x)-p(x)} dx
    probused = 0  # To increase efficiency, we're only to search over a subset of the c values.
    # This will keep track of what probability mass we have covered so far.

    p = math.exp(-eps0)
    expectation = (n-1)*p

    # Now, we are going to iterate over the n/2, n/2-step, n/2+step, n/2-2*steps, ...
    for B in range(1, int(np.ceil(n/step)), 1):
        for s in range(2):
            if s == 0:
                if B==1:
                    upperc = int(np.ceil(expectation+B*step))  # This is stepping up by "step".
                    lowerc = upperc - step
                else:
                    upperc = int(np.ceil(expectation + B * step))  # This is stepping up by "step".
                    lowerc = upperc - step + 1
                if lowerc>n-1:
                    inscope = False
                else:
                    inscope = True
                    upperc = min(upperc, n-1)
            if s == 1:
                lowerc = int(np.ceil(expectation-B*step))
                upperc = lowerc + step - 1
                if upperc<0:
                    inscope = False
                else:
                    inscope = True
                    lowerc = max(0, lowerc)

            if inscope == True:
                cdfinterval = stats.binom.cdf(upperc, n - 1, p) - stats.binom.cdf(lowerc, n - 1, p) + stats.binom.pmf(lowerc, n - 1, p)
            # This is the probability mass in the interval (in Bin(n-1, p))

                if max(deltap, deltaq) > deltaupper:
                    return deltaupper

                if 1 - probused < deltap and 1 - probused < deltaq:
                    if upperbound == True:
                        return max(deltap + 1 - probused, deltaq + 1 - probused)
                    else:
                        return max(deltap, deltaq)

                else:
                    deltap_upperc = onestep(upperc, eps, eps0, True)
                    deltap_lowerc = onestep(lowerc, eps, eps0, True)
                    deltaq_upperc = onestep(upperc, eps, eps0, False)
                    deltaq_lowerc = onestep(lowerc, eps, eps0, False)

                    if upperbound == True:
                        # compute the maximum contribution to delta in the segment.
                        # The max occurs at the end points of the interval due to monotonicity
                        deltapadd = max(deltap_upperc, deltap_lowerc)
                        deltaqadd = max(deltaq_upperc, deltaq_upperc)
                    else:
                        deltapadd = min(deltap_upperc, deltap_lowerc)
                        deltaqadd = min(deltaq_upperc, deltaq_lowerc)

                    deltap = deltap + cdfinterval * deltapadd
                    deltaq = deltaq + cdfinterval * deltaqadd

                probused = probused + cdfinterval  # updates the mass of C covered so far

    return max(deltap, deltaq)


# #if UL=1 then produces upper bound, else produces lower bound.
def numericalanalysis(n, epsorig, delta, num_iterations, step, upperbound):
    '''
    Empirically computes the privacy guarantee of achieved by shuffling n eps0-DP local reports.
    num_iterations = number of steps of binary search, the larger this is, the more accurate the result
    If upperbound=True then this produces an upper bound on the true shuffled eps, and if upperbound=False then it produces a lower bound.
    '''
    # in order to speed things up a bit, we start the search for epsilon off at the theoretical upper bound.
    if epsorig < math.log(n / (16 * math.log(4 / delta))):
        # checks if this is a valid parameter regime for the theoretical analysis.
        # If yes, uses the theoretical upper bound as a starting point for binary search
        epsupper = closedformanalysis(n, epsorig, delta)
    else:
        epsupper = epsorig

    def deltacompinst(eps, delta):
        return deltacomp(n, epsorig, eps, delta, step, upperbound)

    return binarysearch(deltacompinst, delta, num_iterations, epsupper)


# ===========/THEORY/========
def closedformanalysis(n, epsorig, delta):
    '''
    Theoretical computation the privacy guarantee of achieved by shuffling n eps0-DP local reports.
    '''
    if epsorig > math.log(n / (16 * math.log(4 / delta))):
        print("This is not a valid parameter regime for this analysis")
        return epsorig
    else:
        a = 8 * (math.exp(epsorig) * math.log(4 / delta)) ** (1 / 2) / (n) ** (1 / 2)
        c = 8 * math.exp(epsorig) / n
        e = math.log(1 + a + c)
        b = 1 - math.exp(-epsorig)
        d = (1 + math.exp(-epsorig - e))
        return math.log(1 + (b / d) * (a + c))
