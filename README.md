# VANESA
**Cumulative Multiple-Niching Evolution Strategy Algorithm with Adaptive Sampling and BFGS Post-Processing.**
**Version 2.0**

Goal: Obtain LOCAL optima.

This is a second-generation prototype code. It is still based on the old CMNGA code, hence to provide continuity, the repo is still named CMNGA.
I switched from GA to ES as it provided better results for out testing function in fewer steps. There are many changes to the previous code:
0) GA -> ES.
1) BFGS post-processing of the found optima. Initial Hessian is approximated by the covariance matrix.
2) Compartmentalization of the "five peaks function" and inclusion of further functions (incl. Rastrigin's) for further testing.
3) Latin Hypercube Sampling to avoid "clumping". NOTE: There is likely no difference to absolutely random sampling. Or none perceived in our 2D test.
4) Calculation of the Chi critical value. NOTE: The tabulation was good, but the change leaves the code a) easier to expand to higher dimensions, b) easy to switch the percentile, c) I have not found any Fortran code that would include both. The precision is decent enough. Read the comments!
5) License text file.

TO DO:
1) I still have to go over the populations and assignments. The code correctly converges to 4 highest peaks at sufficient samplings (~15 000 function evalfs), at low samplings it either finds the one broadest, or two highest. The 5th peak is drowned in the 'noise'.
2) MPI implementation—split populations.

PREREQUISITES:
1) BLAS/LAPACK = we use single-precision; see module blas_external

IDEAS:
1) As we use SP = GPU acceleration?

*Use at your own risk. Citations to the used literature are added continuously within the code.*
