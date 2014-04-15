PIPM: Perturbed Interior Point Method for LP
===============================================

**Description**:
An infeasible primal-dual path-following interior point method with 
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
        options - struct, optional input, controlling all parameters
        	.maxIter        Maximum number of iterations allowed
        	.tol            Convergence tolerance
        	.mu_cap         Threshold value for mu
        	.cutoff         Threshold value for cutoff
        	.iPer           Initial perturbations

        	.actvPredStrtgy String, determine the strategy of active-set prediction
                     		By default it's 'conservCutoff'.
                     		'simple' - simple cutoff
                     		'conservCutoff' - three-set with cutoff
                     		'conservIdFunc' - three-set with idFunc
                     		'conservIndica' - three-set with indicator

        	.doCrossOver 	Controls whether or not perform crossover to
                    		simplex after ipm iterations.
                    		Default value 1.
                          	0 - No
                          	1 - Yes

        	.verbose = 2    Controls how much information to display.
                 	        0   - Nothing
                          	1   - Only optimal information
                          	2   - every iterations
                          	>=3 - All information

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
