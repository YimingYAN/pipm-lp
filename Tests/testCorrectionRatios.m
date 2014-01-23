% FUNCTION testCorrectionRatios
%
% -------------------------------------------------------------------------
%
% This function is used to compare the correction ratios for the algorithms
% with and without perturbations.
%
% -------------------------------------------------------------------------
%
% The results generated by this scriptis are presented in our paper
%
% "Active-set prediction for interior point methods using perturbations"
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
% Notes:
%   1. In the paper, we denote the perturbed algorithm as Algorithm 6.1 and
%   the unperturbed algorithm as Algorithm 6.4.
%
% -------------------------------------------------------------------------
%
% 22 Jan 2014
% Yiming Yan
% University of Edinburgh

%% ----------------- Main Func -------------------- %%
function testCorrectionRatios
clear all;
close all;
clc;

%% Setup
[Type, numTestProb, params_per, params_unper] = setup_correctionRatio;

% -------------------------------------------------------------------------
stopAtRangeL = 8;
stopAtRangeU = 12;

% Options for plots
% Legends = {'With perturbations' 'Without perturbations'};
Legends = {'Algorithm 6.1' 'Algorithm 6.4'};
fileName = ['correction_ratio_test_' Type];

%% Run the test
fprintf('\n2. Start the %s test...\n', Type);
fprintf('\n================================= Correction Ratio Tests =================================\n');
printHeader;

counter = 1;      % Counter for outter loop

falsePrediction  = zeros(stopAtRangeU - stopAtRangeL, 2);
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
    
    fpr1 = 0; mpr1 = 0; cr1 = 0;
    fpr2 = 0; mpr2 = 0; cr2 = 0;
    
    avgRes1 = 0; avgRes2 = 0;
    
    Avgm = 0; Avgn = 0;
    
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
        
        [m, n] = size(A);
        %% Solve the problem using linprog (ipm) and get the actual actv
        [actualActv, exitflag] = solveLinprog(A,b,c);
        
        % Skip the test problem if linprog does not converge
        if exitflag ~= 1
            skipped = skipped + 1;
            continue;
        end
        
        %% Solve the problem using pipm-lp (Perturbed alg)
        params_per.maxIter = k;
        per = pipm(A,b,c,params_per); per.solve;
        
        params_unper.maxIter = k;
        unper = pipm(A,b,c,params_unper); unper.solve;
        
        %% Accumulate the correction ratios
        [tpfpr1, tpmpr1, tpcr1] = ...
            getCorrectionRatio(per.getActv, actualActv);
        [tpfpr2, tpmpr2, tpcr2] = ...
            getCorrectionRatio(unper.getActv, actualActv);
        
        fpr1 = fpr1+tpfpr1; mpr1 = mpr1+tpmpr1; cr1 = cr1+tpcr1;
        fpr2 = fpr2+tpfpr2; mpr2 = mpr2+tpmpr2; cr2 = cr2+tpcr2;
        
        avgRes1 = per.getIPMResidual + avgRes1;
        avgRes2 = unper.getIPMResidual   + avgRes2;
        
        Avgm = Avgm+m; Avgn = Avgn+n;
        
        %% Increment the counter
        i = i+1;
    end
    
    %% Get the averages
    falsePrediction(counter,:)  = [fpr1 fpr2]./numTestProb;
    missedPrediction(counter,:) = [mpr1 mpr2]./numTestProb;
    correctionR(counter,:)      = [cr1 cr2]./numTestProb;
    
    avgResidual(counter,:)      = [avgRes1 avgRes2]./numTestProb;
    
    Avgm = round(Avgm/numTestProb); Avgn = round(Avgn/numTestProb);
    
    printContent(k, counter, Avgm, Avgn, falsePrediction,...
        missedPrediction, correctionR, avgResidual);
    
    counter = counter+1;
end

fprintf('\n\tTotal number of probs skipped: %d\n', skipped);

clearvars A b c lb ub m n i k xsol exitflag...
    Prob per unper cr* fpr* mpr* tp*...
    avgRes1 avgRes2 actualActv counter;
save( ['correction_ratio_test_' Type '.mat'] );

%% Output the result
fprintf('\n3. Output the result...\n');
range = stopAtRangeL : stopAtRangeU;

plotCorrectionRatios(falsePrediction, missedPrediction,...
    correctionR, avgResidual, range, Legends,fileName);

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

numTestProb  = 100; % set to 10 for demo. 100 for real test.

% With perturbations
parameters_per.verbose          = 0;
parameters_per.iPer             = 1e-02;
parameters_per.actvPredStrtgy   = actvPredStrtgy;
parameters_per.doCrossOver      = 0;

% Without perturbations
parameters_unper.verbose        = 0;
parameters_unper.iPer           = 0;
parameters_unper.actvPredStrtgy = actvPredStrtgy;
parameters_unper.doCrossOver    = 0;
end

%% Function used to solve the LP using linprog
function [actualActv, exitflag] = solveLinprog(A,b,c)

n  = size(A, 2);   A  = full(A);
lb = zeros(n, 1);  ub = inf*ones(n, 1);

% options = optimset('LargeScale', 'off', 'Algorithm','simplex','Display','off');
% options = optimset('Algorithm','active-set','Display','off');
options = optimset('Algorithm', 'interior-point','Display','off');

[xsol,~,exitflag] = linprog(c,[],[],A,b,lb,ub,[],options);

% Get actual actv
actualActv = find(xsol<1e-05);
end

%% Function used to calculate the correction ratios
function [falsePrediction, missedPrediction, cr] = ...
    getCorrectionRatio(predictedActv,actualActv)

% actualActv is empty?
if isempty(actualActv)
    if ~isempty(predictedActv)
        falsePrediction = 1;
    else
        falsePrediction = 0;
    end
    
    missedPrediction = 0;
    cr = 1-falsePrediction;
else
    M = length(union(predictedActv,actualActv));
    
    falsePrediction =...
        length(setdiff(predictedActv,actualActv))/M;
    missedPrediction = ...
        length(setdiff(actualActv,predictedActv))/M;
    
    cr = 1-falsePrediction - missedPrediction;
end


end

%% This function is used to plot the correction ratios
function plotCorrectionRatios(falsePrediction, missedPrediction,...
    correctionR, avgResidual, range, Legends, fileName)
clf;
lineStyle = {'-ro' '--b*'};

h(1) = subplot(2,2,1);
plot(range,falsePrediction(:,1),lineStyle{1},...
    range,falsePrediction(:,2),lineStyle{2});
title('False Prediction Ratio');

h(2) = subplot(2,2,2);
plot(range,missedPrediction(:,1),lineStyle{1},...
    range,missedPrediction(:,2),lineStyle{2});
title('Missed Prediction Ratio');

h(3) = subplot(2,2,3);
plot(range,correctionR(:,1),lineStyle{1},...
    range,correctionR(:,2),lineStyle{2});
title('Correction Ratio');

h(4) = subplot(2,2,4);
plot(range,log10(avgResidual(:,1)),lineStyle{1},...
    range,log10(avgResidual(:,2)),lineStyle{2})
title('Average Relative Residual (log10)')

% linkaxes(h,'xy');
set(h,'XTick',range);
set(h,'XLim',[range(1) range(end)+0.1]);
set(h(1:3),'YLim',[-0.1 1.1]);
set(h(4), 'YLim',[floor(min(min(log10(avgResidual)))) ceil(max(max(log10(avgResidual))))]);
set(h,'XGrid','on','YGrid','on');
set(h(4),'YTick', floor(min(min(log10(avgResidual)))):1:ceil(max(max(log10(avgResidual)))));
hleg = legend(Legends,'Orientation','horizontal');
p =get(hleg,'Position');
p(1) = 0.5-0.5*p(3); p(2)= 0.02;
set(hleg,'Position',p);
print('-depsc',fileName);
[result,msg] = eps2pdf([fileName '.eps']);

end

%% Print iterative info
function printHeader
fprintf('\n%4s | %11s | %7s %7s %7s %9s | %7s %7s %7s %9s\n',...
    'Iter', '[m,n]  ', 'FPR_PER', 'MPR_PER', 'CR_PER', 'RES_PER',...
    'FPR_UNP', 'MPR_UNP', 'CR_UNP', 'RES_UNP');
end

function printContent(k, counter, Avgm, Avgn, falsePrediction,...
    missedPrediction, correctionR, avgResidual)
fprintf('%4d | [%4d %4d] | %7.2f %7.2f %7.2f %9.2e | %7.2f %7.2f %7.2f %9.2e\n',...
    k,Avgm,Avgn,...
    falsePrediction(counter, 1), missedPrediction(counter, 1), ...
    correctionR(counter, 1),      avgResidual(counter, 1), ...
    falsePrediction(counter, 2), missedPrediction(counter, 2), ...
    correctionR(counter, 2),      avgResidual(counter, 2) );
end

