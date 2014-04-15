Perturbed Interior Point Method (PIPM)
===============================================

An infeasible primal-dual path-following interior point method with 
controlled perturbations for Linear Programming.

The main solver PIPM (/src/pipm.m) only accepts LP problems 
in the standard form, i,e,
```
        min c'*x s.t. Ax = b, x>=0.        
```

The future version may also accept QP problems, but currently it has not been implemented.


**Syntax**:
```
        p = pipm( A, b, c);             p.solve;
        p = pipm( A, b, c, options);    p.solve;
```

**Input**:

```
        (A,b,c) - problem data
        options - struct, optional input, controlling all parameters
        	.maxIter        Maximum number of iterations allowed
                            Default value 100
        	                
        	.tol            Convergence tolerance
                            Default value 1e-06
        	                
        	.mu_cap         Threshold value for mu
        	                Default value 1e-03
        	                
        	.cutoff         Threshold value for cutoff
        	                Default value 1e-05
        	                
        	.iPer           Initial perturbations
                            Default value 1e-02
                                
        	.actvPredStrtgy String, determine the strategy of active-set prediction
                            Default value 'conservCutoff'.
                     		'simple' - simple cutoff
                     		'conservCutoff' - conservative strategy with cutoff
                     		'conservIdFunc' - conservative strategy with idFunc
                     		'conservIndica' - conservative strategy with indicators

        	.doCrossOver 	Controls whether or not perform crossover to
                    		simplex after ipm iterations.
                    		Default value 1
                          	0 - No
                          	1 - Yes

        	.verbose        Controls how much information to display.
                            Default value 2
                 	        0   - Nothing
                          	1   - Only optimal information
                          	2   - Iterative information
                          	>=3 - All information

```

**PIPM as normal infeasible primal-dual path-following IPM**:
```
        options.iPer = 0;               % No perturbations
        options.doCrossOver = 0;        % Disable crossover to simplex
        options.mu_cap = 1e-09;         % Push duality gap to sufficiently small
        
        p = pipm( A, b, c, options);
        p.solve;
```

**Example**: 

Please see examples.m in folder Examples


**Author**:  

Yiming Yan @ University of Edinburgh


**Reference**:

Coralia Cartis and Yiming Yan, 
[Active-set prediction for interior point methods using controlled perturbations](http://www.maths.ed.ac.uk/~yan/research/papers/aiipm.pdf),
submitted, April 2014

- - -

Have Fun! :)

Yiming Yan,
University of Edinburgh

April 15, 2014
