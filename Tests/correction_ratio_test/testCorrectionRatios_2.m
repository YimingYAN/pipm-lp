function testCorrectionRatios_2
% FUNCTION testCorrectionRatios_2
% This function is programmed to verfy Theorem ?.? in our paper [1].
%
% -------------------------------------------------------------------------
%
% Plot the false-prediction, missed-prediction and correction ratios for
% for the following 4 cases:
%
% 1. Predicted original active-set from perturbed alg  \
%                           V.S.                       | Alg 6.1 - Splx
%    Actual original activ-set from linprog (simplex)  /
%
% 2. Predicted original active-set from perturbed alg  \
%                           V.S.                       | Alg 6.1 - IPM
%    Actual original activ-set from linprog (ipm)      /
%
% 3. Predicted perturbed active-set from perturbed alg \
%                           V.S.                       | Alg 6.2_P - Splx_P
%    Actual perturbed activ-set from linprog (simplex) /
%
% 4. Predicted perturbed active-set from perturbed alg \
%                           V.S.                       | Alg 6.2_P - IPM_P
%    Actual perturbed activ-set from linprog (ipm)     /
%
% -------------------------------------------------------------------------
%
% How to use this script:
%   1. Run the script
%   2. Choose the test set
%   3. Check the results (.mat, .pdf, .eps)
%
% -------------------------------------------------------------------------
%
% [1] Active-set prediction for interior point methods using perturbations
%
%
% Date        :  22 Feb 2014
% Author      :  Yiming Yan
% Affiliation :  University of Edinburgh

close all;
clc;

%% Setup
[Type, numTestProb, params_oac_per, params_pac_per] = setup_correctionRatio;

% -------------------------------------------------------------------------
stopAtRangeL = 6;
stopAtRangeU = 18;

% Options for plots
Legends = { 'Alg 6.1 - Splx'  'Alg 6.1 - IPM'  'Alg 6.2\_P - Splx\_P' 'Alg 6.2\_P - IPM\_P' };
% lineStyle = {'-ro' '--b*' '-kv' '--m+'};
colors = {'r' 'b' 'y' 'g'};
lineStyles = {'-' '--' '-' '--'};
markers = {'o' '*' 'v' '+'};
fileName = ['correction_ratio_test_2_' Type];

%% Run the test
fprintf('\n2. Start the %s test...\n', Type);
fprintf('\n================================= Correction Ratio Tests =================================\n');
printHeader;

counter = 1;      % Counter for outter loop

falsePrediction  = zeros(stopAtRangeU - stopAtRangeL, 4);
missedPrediction = falsePrediction;
correctionR      = falsePrediction;
avgResidual      = falsePrediction;

skipped = 0;

for k=stopAtRangeL:1:stopAtRangeU
    %% Reset the random number generator
    rng('default');
    rng(1);
    
    % Initialize data
    i=1;
    
    fpr_per_splx = zeros(numTestProb,1);
    mpr_per_splx = fpr_per_splx; cr_per_splx  = fpr_per_splx;
    
    fpr_unp_splx = fpr_per_splx;
    mpr_unp_splx = fpr_per_splx; cor_unp_splx = fpr_per_splx;
    
    fpr_per_ipm  = fpr_per_splx;
    mpr_per_ipm  = fpr_per_splx; cor_per_ipm  = fpr_per_splx;
    
    fpr_unp_ipm  = fpr_per_splx;
    mpr_unp_ipm  = fpr_per_splx; cor_unp_ipm  = fpr_per_splx;
    
    avgRes_per = fpr_per_splx; avgRes_unper = fpr_per_splx;
    
    Avgm = zeros(numTestProb,1); Avgn = zeros(numTestProb,1);
    
    while i<=numTestProb
        %% Generate random problem
        switch lower(Type)
            case 'random'
                [A, b, c] = ...
                    generateRandomProb('m_min',10,'m_max',200,...
                    'n_min',20,'n_max',500);
            case 'random_degen'
                [A, b, c] =...
                    generateDegenProb('m_min',10,'m_max',200,...
                    'n_min',20,'n_max',500);
            otherwise
                return;
        end
        
        
        %% Solve for 'actual active-set' of the original problem
        [m, n] = size(A);
        
        %  Get the actual original actv from linprog (simplex) 
        [ actualActv_splx,   exitflag_per ] = solveLinprog(A, b, c,'splx');
        
        %  Get the actual original actv from linprog (ipm)
        [ actualActv_ipm, exitflag_unper  ] = solveLinprog(A, b, c,'ipm' );
        
        % Skip the test problem if linprog does not converge
        if exitflag_per ~= 1 || exitflag_unper ~= 1
            skipped = skipped + 1;
            continue;
        end
        
        %% Predict the original actv using pipm-lp (Perturbed alg)
        params_oac_per.maxIter = k;
        oac_per = pipm(A,b,c,params_oac_per); oac_per.solve;
        
        %% Build the perturbed problem
        bp = b + A * oac_per.perturbations.lambda;
        cp = c + oac_per.perturbations.lambda;
        
        %% Predict the perturbed actv using pipm-lp (solve the euqavilent perturbed problems with iPer = 0)
        params_pac_per.maxIter = k;
        pac_per = pipm(A, bp, cp, params_pac_per); pac_per.solve;
        
        %% Find the actual perturbed actv using
        [actualActv_per_splx,   ~  ] = solveLinprog(A, bp, cp, 'splx');
        [actualActv_per_ipm,    ~  ] = solveLinprog(A, bp, cp, 'ipm' );
        
        %% Get the correction ratios
        
        % Predicted original actv from pipm-lp VS original actv from simplex
        [fpr_per_splx(i), mpr_per_splx(i), cr_per_splx(i)] = ...
            getCorrectionRatio(oac_per.getActv, actualActv_splx);
        
        % Predicted original actv from pipm-lp VS original actv from ipm
        [fpr_per_ipm(i), mpr_per_ipm(i), cor_per_ipm(i)] = ...
            getCorrectionRatio(oac_per.getActv, actualActv_ipm);
        
        % Predicted perturbed actv from pipm-lp VS original perturbed actv
        % from simplex
        [fpr_unp_splx(i), mpr_unp_splx(i), cor_unp_splx(i)] = ...
            getCorrectionRatio(pac_per.getActv, actualActv_per_splx);
        
        % Predicted perturbed actv from pipm-lp VS original perturbed actv
        % from ipm
        [fpr_unp_ipm(i), mpr_unp_ipm(i), cor_unp_ipm(i)] = ...
            getCorrectionRatio(pac_per.getActv, actualActv_per_ipm);
        
        avgRes_per(i)   = oac_per.getIPMResidual;
        avgRes_unper(i) = pac_per.getIPMResidual;
        
        Avgm(i) = m; Avgn(i) = n;
        
        %% Increment the counter
        i = i+1;
    end
    
    %% Get the averages
    falsePrediction(counter,:)  = mean( [ fpr_per_splx  fpr_per_ipm  fpr_unp_splx fpr_unp_ipm ] );
    missedPrediction(counter,:) = mean( [ mpr_per_splx  mpr_per_ipm  mpr_unp_splx mpr_unp_ipm ] );
    correctionR(counter,:)      = mean( [ cr_per_splx   cor_per_ipm  cor_unp_splx cor_unp_ipm ] );
    
    avgResidual(counter,:)      = mean([avgRes_per avgRes_per avgRes_unper avgRes_unper]);
    
    Avgm = round(mean(Avgm));   Avgn = round(mean(Avgn));
    
    printContent(k, counter, Avgm, Avgn, falsePrediction,...
        missedPrediction, correctionR, avgResidual);
    
    save
    
    counter = counter+1;
end

fprintf('\n\tTotal number of probs skipped: %d\n', skipped);

clearvars A b c lb ub m n i k xsol exitflag...
    Prob per unper cr* fpr* mpr* tp*...
    avgRes1 avgRes2 actualActv counter;

%% Output the result
fprintf('\n3. Output the result...\n');
range = stopAtRangeL : stopAtRangeU;

save( ['correction_ratio_test_2_' Type '.mat'] );

plotCorrectionRatios_4lines(falsePrediction, missedPrediction,...
    correctionR, avgResidual, range, Legends, fileName,...
    colors, lineStyles, markers);

fprintf('DONE.\n');
fprintf('Pls check the file %s for the plot.\n', [fileName '.pdf']);
end

%% ----------------- Main Func End ----------------- %%



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

numTestProb  = 5; % set to 10 for demo. 100 for real test.

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

%% Function used to solve the LP using linprog
function [actualActv, exitflag] = solveLinprog(A,b,c, alg)

n  = size(A, 2);   A  = full(A);
lb = zeros(n, 1);  ub = inf*ones(n, 1);

switch lower(alg)
    case 'splx'
        options = optimset('LargeScale', 'off', 'Algorithm','simplex','Display','off');
    case 'ipm'
        options = optimset('Algorithm', 'interior-point','Display','off');
    case 'actv-set'
        options = optimset('Algorithm','active-set','Display','off');
    otherwise
        error('solveLinprog: no such a solver');
end

[xsol,~,exitflag] = linprog(c,[],[],A,b,lb,ub,[],options);

% Get actual actv
actualActv = find(xsol<1e-05);
end

%% Print iterative info
function printHeader
fprintf('\n%4s | %11s | %7s %7s %7s %9s | %7s %7s %7s %9s | %7s %7s %7s %9s | %7s %7s %7s %9s\n',...
    'Iter', '[m,n]  ',...
    'F_oA_PS', 'M_oA_PS', 'C_oA_PS', 'R_oA_PS',...
    'F_oA_PI', 'M_oA_PI', 'C_oA_PI', 'R_oA_PI',...
    'F_pA_PS', 'M_pA_PS', 'C_pA_PS', 'R_pA_PS',...
    'F_pA_PI', 'M_pA_PI', 'C_pA_PI', 'R_pA_PI');
end

function printContent(k, counter, Avgm, Avgn, falsePrediction,...
    missedPrediction, correctionR, avgResidual)
fprintf('%4d | [%4d %4d] | %7.2f %7.2f %7.2f %9.2e | %7.2f %7.2f %7.2f %9.2e | %7.2f %7.2f %7.2f %9.2e | %7.2f %7.2f %7.2f %9.2e\n',...
    k,Avgm,Avgn,...
    falsePrediction(counter, 1), missedPrediction(counter, 1), ...
    correctionR(counter, 1),      avgResidual(counter, 1), ...
    falsePrediction(counter, 2), missedPrediction(counter, 2), ...
    correctionR(counter, 2),      avgResidual(counter, 2),...
    falsePrediction(counter, 3), missedPrediction(counter, 3), ...
    correctionR(counter, 3),      avgResidual(counter, 3),...
    falsePrediction(counter, 4), missedPrediction(counter, 4), ...
    correctionR(counter, 4),      avgResidual(counter, 4) );

end

