% This script shows several examples, demonstrating how to use solver PIPM.
clear;
clc;

%% Simple LP 
A = [1 1 -1  0;
     1 2  0 -1];
b = [ 3; 4];
c = [ 2; 1 ; 0; 0];

p = pipm(A,b,c);
p.solve;

clear;
%% Example on Randomly Generated Primal Nondegenrate LP
[A, b, c] = generateRandomProb;
p = pipm(A,b,c);
p.solve;

clear
%% Example on Randomly Generated Primal Dengenerate LP
[A, b, c] = generateDegenProb;
p = pipm(A,b,c);
p.solve;

clear
%% Example on changing initial perturbations
A = [1 1 -1  0;
     1 2  0 -1];
b = [ 3; 4];
c = [ 2; 1 ; 0; 0];

parameters_input.iPer = 1e-03;

p = pipm(A,b,c,parameters_input);
p.solve;

clear;
%% Example on changing mu_cap
[A, b, c] = generateRandomProb;

parameters_input.mu_cap = 1e-02;

p = pipm(A,b,c,parameters_input);
p.solve;

clear;

%% Example on changing the nuber of maximum ipm iterations allowed,  maxIter
[A, b, c] = generateRandomProb;

parameters_input.maxIter = 2;

p = pipm(A,b,c,parameters_input);
p.solve;

clear;

%% Example on chaning verbose
[A, b, c] = generateRandomProb;

parameters_input.verbose = 0;
p = pipm(A,b,c,parameters_input);
p.solve; 

parameters_input.verbose = 1;
p = pipm(A,b,c,parameters_input);
p.solve; 

parameters_input.verbose = 2;
p = pipm(A,b,c,parameters_input);
p.solve; 

parameters_input.verbose = 3;
p = pipm(A,b,c,parameters_input);
p.solve; 

clear;

%% Exmaple on choosing different actvPredStrategies
[A, b, c] = generateRandomProb;
parameters_input.mu_cap = 1e-03;


% simple
parameters_input.actvPredStrtgy = 'simple';
p = pipm(A,b,c,parameters_input);
p.solve;

% conservCutoff
parameters_input.actvPredStrtgy = 'conservCutoff';
p = pipm(A,b,c,parameters_input);
p.solve;

% conservIdFunc
parameters_input.actvPredStrtgy = 'conservIdFunc';
p = pipm(A,b,c,parameters_input);
p.solve;

% wrong setup
parameters_input.actvPredStrtgy = 'wrongSetup';
p = pipm(A,b,c,parameters_input);
try 
     p.solve;
catch err
     fprintf('%s\n',err.message)
end

clear;

%% chekc if Netlib problems work
load AFIRO;
parameters_input.mu_cap = 1e-03;
p = pipm(A,b,c,parameters_input);
p.setProbName('AFIRO');

p.solve

clear;
