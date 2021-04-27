# For licensing see accompanying LICENSE file.
# Copyright (C) 2020 Apple Inc. All Rights Reserved.

import computeamplification as CA

#number of local reports
n=10**5
#local epsilon
epsorig = 4
#desired delta of final shuffled privacy guarantee
delta = 10**(-6)
#number of iterations of binary search. The higher T is, the more accurate the result
num_iterations = 10
#This is a parameter of the empirical analysis computation that can be tuned for efficiency. The larger step is, the less accurate the result, but more efficient the algorithm.
step = 100

#There are 2 main functions, empiricalanalysis and theoryanalysis.
#empiricalanalysis computes the shuffled privacy guarantee empirically. The 1 and 0 correspond to returning either an upper bound on the privacy guarantee, or a lower bound.
empirical_upperbound = CA.empiricalanalysis(n, epsorig, delta, num_iterations, step, True)
empirical_lowerbound = CA.empiricalanalysis(n, epsorig, delta, num_iterations, step, False)

print("Shuffling", n, epsorig, "-DP local randomizers results is (eps, ", delta, ")-DP in the shuffle model for eps between", empirical_lowerbound, "and", empirical_upperbound)

#theoryanalysis computes the privacy amplification based on our theoretical analysis.
theoretical = CA.theoryanalysis(n, epsorig, delta)

print("According to our theoretical analysis, shuffling", n, epsorig, "-DP local randomizers results is (eps, ", delta, ")-DP in the shuffle model where eps is at most ", theoretical)
