# Simulation codes for testing validation of PGA's global optimality

Please directly run mainfile.m to obtain the test results.

Qp_gen: function for generating data Q and p.

PGAforQP: function for iteratively solving a convex quadratic programming model without the m-sparse constraint by PGA.

L0QuadProg: function for exhaustively enumerating support sets to obtain the globally optimal solution of model (3.7).

prox_mpsparse: function for computing the proximity operator of the indicator function of m-sparse constraint.

PGAforL0QuadProg: function for iteratively solving model (3.7) in the paper by PGA.
