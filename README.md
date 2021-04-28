## Hiding Among the Clones: Privacy Amplification by Shuffling

This software project accompanies the research paper, Hiding Among the Clones: A Simple and Nearly Optimal Analysis of Privacy Amplification by Shuffling (https://arxiv.org/abs/2012.12803).

Recent work of Erlingsson, Feldman, Mironov, Raghunathan, Talwar, and Thakurta [EFMRTT19] demonstrates that random shuffling amplifies differential privacy guarantees of locally randomized data. Such amplification implies substantially stronger privacy guarantees for systems in which data is contributed anonymously [BEMMRLRKTS17] and has lead to significant interest in the shuffle model of privacy [CSUZZ19,EFMRTT19].  We show that random shuffling of n data records that are input to ε_0-differentially private local randomizers results in an (O((1−e^(ε_0))sqrt(e^(ε_0)log(1/δ)/n)),δ)-differentially private algorithm. This significantly improves over previous work and achieves the asymptotically optimal dependence in ε_0. Our result is based on a new approach that is simpler than previous work and extends to approximate differential privacy with nearly the same guarantees. Our work also yields an empirical method to derive tighter bounds the resulting ε and we show that it gets to within a small constant factor of the optimal bound. As a direct corollary of our analysis, we derive a simple and asymptotically optimal algorithm for discrete distribution estimation in the shuffle model of privacy. We also observe that our result implies the first asymptotically optimal privacy analysis of noisy stochastic gradient descent that applies to sampling without replacement.

This repository contains a python implementation of our main amplification bound from shuffling purely differentially private local randomizers. 

## Documentation

Appendix E of https://arxiv.org/abs/2012.12803 contains details of our implementation. Some function names differ, "B" in appendix is "onestep" in code and "M_1" in appendix is "deltacomp" in code.

## Getting Started 

example.py contains an example showing how to compute amplification bounds with this code. Given the parameters set as default in example.py, the output of this script should be the following two lines:

Shuffling 100000 4 -DP local randomizers results is (eps,  1e-06 )-DP in the shuffle model for eps between 0.1675385583317841 and 0.172790550755978

According to our closed form analysis, shuffling 100000 4 -DP local randomizers results is (eps,  1e-06 )-DP in the shuffle model where eps is at most  0.5378040242374512