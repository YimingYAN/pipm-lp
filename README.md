Perturbed Interior Point Method for LP (PIPM)
===============================================

An infeasible primal-dual path-following interior point method with 
controlled perturbations.

The main solver PIPM (/src/pipm.m) only accepts LP problems 
in the standard form, i,e,
```
        min c'*x s.t. Ax = b, x>=0.        
```

The future version may also accept QP problems, but currently it has not been implemented.

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
                     		Default value 'conservCutoff'.
                     		'simple' - simple cutoff
                     		'conservCutoff' - three-set with cutoff
                     		'conservIdFunc' - three-set with idFunc
                     		'conservIndica' - three-set with indicator

        	.doCrossOver 	Controls whether or not perform crossover to
                    		simplex after ipm iterations.
                    		Default value 1
                          	0 - No
                          	1 - Yes

        	.verbose        Controls how much information to display.
                            Default value 2
                 	        0   - Nothing
                          	1   - Only optimal information
                          	2   - every iterations
                          	>=3 - All information

```
---------------------------------------------------------------------------
**PIPM as normal infeasible primal-dual path-following IPM**:
```
        options.iPer = 0;               % No perturbations
        options.doCrossOver = 0;        % Disable crossover to simplex
        options.mu_cap = 1e-09;         % Push duality gap to sufficiently small
        
        p = pipm( A, b, c, options);
        p.solve;
```
---------------------------------------------------------------------------
**Example**: 

Please see examples.m in folder Examples

---------------------------------------------------------------------------
**Author**:  

Yiming Yan @ University of Edinburgh

---------------------------------------------------------------------------
**Reference**:

Coralia Cartis and Yiming Yan, 
[Active-set prediction for interior point methods using controlled perturbations](http://www.maths.ed.ac.uk/~yan/research/papers/aiipm.pdf),
submitted, April 2014

---------------------------------------------------------------------------

Have Fun! :)

Yiming Yan

University of Edinburgh

April 15, 2014
