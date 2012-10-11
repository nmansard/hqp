hqp
===

Hierachical Quadratic Problem solver

This software implements a MATLAB-based hierarchical least-square quadratic
solver (HQP) based on the paper 

The software is organize as follow:

1. in the main directory, the main functions of the solver:

* hcod computes the HCOD from the problem specification A,b$ and an initial
  guess of the active set.  The function returns h and Y that are manipulated
  by the other function during the active-search loop.

* eHQP_primal computes the primal optimum of the eHQP corresponding to the
  current active set. The function returns the optimum x and the corresponding
  expression in the Y basis y. This corresponds to Alg. 1 of the paper

* step_length corresponds to Alg. 5 of the paper. The function returns the
  step length, along with a Boolean specifying if one constraint needs to be
  activated and the constraint reference.

* eHQP_dual compute the multiplier of one level of the hierarchy,
  as proposed in Alg. 2 of the paper.

* check_mult checks if all the multipliers are positive for lexicographic
  order. The function returns the Boolean ``one multiplier is negative'' along
  with the lowest multiplier reference.

* active_search is the main function: it computes the HCOD and then run the
  search loop. It returns the primal and dual optima and the HCOD.

2. in utils, the secondary functions used by the solver.

3. in utest, some unitary tests that validates the main functionnalities of the
  solver.

4. in utest/utils, a set of secondary functions used by the tests, with no real
  arrangement nor documentation.

All the code should be run from the main directory. Add utest in your path and
launch the unitary tests:

    addpath('utest')
    tsearch

Author
======
Nicolas Mansard, LAAS/CNRS, Toulouse, France

Reference 
=========

* "Hiearchical Quadratic Programming", by Adrien Escande (JRL-Japan/CNRS,
Tsukuba, Japan), Nicolas Mansard (LAAS/CNRS, Toulouse, France) and Pierre-Brice
Wieber (INRIA, Grenoble, France), on-going submission to IJRR.