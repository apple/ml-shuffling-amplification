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
numerical_upperbound = CA.numericalanalysis(n, epsorig, delta, num_iterations, step, True)
numerical_lowerbound = CA.numericalanalysis(n, epsorig, delta, num_iterations, step, False)

print("Shuffling", n, epsorig, "-DP local randomizers results is (eps, ", delta, ")-DP in the shuffle model for eps between", numerical_lowerbound, "and", numerical_upperbound)

#theoryanalysis computes the privacy amplification based on our theoretical analysis.
closedform = CA.closedformanalysis(n, epsorig, delta)

print("According to our closed form analysis, shuffling", n, epsorig, "-DP local randomizers results is (eps, ", delta, ")-DP in the shuffle model where eps is at most ", closedform)
