%% Function needed to setup the test
function [Type, numTestProb, parameters_per, parameters_unper] = setup_correctionRatio()
% Determine which set of problems to test on.
% Choose from the following two values:
% random, random_degen
fprintf('1. Pls choose the test set [1-2]: \n');
fprintf('\t [1]. Random test (primal nondegenerate)\n');
fprintf('\t [2]. Random test (primal-dual degenerate)\n');
usrinput_type = input('Your choice here [1-2]: ');

if usrinput_type == 1
    Type = 'random';
elseif usrinput_type == 2
    Type = 'random_degen';
else
    error('testCorrectionRatios: please choose a number from the above list');
end

% Determine which active-set prediction strategy to use.
% In the paper, we mainly show the results of using a constant as threshold
% value ('conservCutoff'). In the last part the paper, we also mentioned
% the use of identification fucntion ('conservIdFunc'). Both strategies have
% been implemented.
actvPredStrtgy = 'conservCutoff'; % Default value conservCutoff
% Alternative: conservIdFunc

numTestProb  = 100; % set to 10 for demo. 100 for real test.

% With perturbations
parameters_per.verbose          = 0;
parameters_per.iPer             = 1e-02;
parameters_per.actvPredStrtgy   = actvPredStrtgy;
parameters_per.doCrossOver      = 0;
parameters_per.mu_cap           = 1e-09;     % terminate by mu_cap and tol
parameters_per.tol              = 1e-09;     % to aviod ill-conditioning

% Without perturbations
parameters_unper.verbose        = 0;
parameters_unper.iPer           = 0;
parameters_unper.actvPredStrtgy = actvPredStrtgy;
parameters_unper.doCrossOver    = 0;
parameters_unper.mu_cap         = 1e-09;     % terminate by mu_cap and tol
parameters_unper.tol            = 1e-09;     % to aviod ill-conditioning
end
