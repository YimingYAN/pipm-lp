clc;
clear

tol = 1e-010

rng('default');
rng(1);

for i =1:10
    
    
    [A, b, c] = ...
        generateRandomProb('m_min',10,'m_max',200,...
        'n_min',20,'n_max',500);
    
    %% With perturbations
    params_per.verbose        = 0;
    params_per.iPer           = 1e-02;
    params_per.actvPredStrtgy = 'conservCutoff';
    params_per.doCrossOver    = 0;
    params_per.mu_cap         = tol;    % Avoid being terminated by mu_cap
    params_per.tol            = tol;
    
    
    %% Without perturbations
    params_unper.verbose        = 0;
    params_unper.iPer           = 0;
    params_unper.actvPredStrtgy = 'conservCutoff';
    params_unper.doCrossOver    = 0;
    params_unper.mu_cap         = tol;
    params_unper.tol            = tol;
    
    
    %% Per
    disp('per ---')
    per = pipm(A,b,c,params_per); per.solve; 
    disp('       --- per')
    %% Unper
    disp('unper ---')
    unper = pipm(A,b,c,params_unper); unper.solve;
    disp('       --- unper')
end