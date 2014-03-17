PIPM: Perturbed Interior Point Methods for LP
===============================================

**Description**:
An infeasible primal-dual path-following inteiror point method with 
controlled perturbations.

The main solver PIPM (/src/pipm.m) only accepts LP problems 
in the standard form, i,e,
```
        min c'*x s.t. Ax = b, x>=0.        
```
The future version may also accept QP problems, but currently it has not 
been implemented.
---------------------------------------------------------------------------
**Syntax**:
```
        p = pipm( A, b, c, options);
        p.solve;
```
---------------------------------------------------------------------------
**Input**: 
```
        (A,b,c) - problem data
        options: optional input, controlling all parameters
```
---------------------------------------------------------------------------
**Example**:
        Please see examples.m in folder Examples
---------------------------------------------------------------------------
**Author**: Yiming Yan @ University of Edinburgh
---------------------------------------------------------------------------
Have Fun! :)

Yiming Yan
University of Edinburgh
September 20, 2013
