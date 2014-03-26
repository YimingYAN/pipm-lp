function testCrossOver2Simplex_IF_ID_CF(Type)
% FUNCTION testCrossOver2Simplex_IF_ID_CF
%
% -------------------------------------------------------------------------
%
% This function is used to compare simplex iteration counts for the
% infeasible pathfollowing ipm with identification function (IF),
% indicators (ID) and cutoff (CF).
% -------------------------------------------------------------------------
%
% Input:
%       Type: random, random_degen, netlib
%
% How to use this script:
%   1. Run the script
%   2. Choose the test set
%   3. Check the results (.mat, .pdf, .eps)
%
% -------------------------------------------------------------------------
% 24 March 2013
% Yiming Yan
% University of Edinburgh

%% %%%%% %%%%%%% %%%%%%% --- Main Func --- %%%%%%% %%%%%%% %%%%% %%
close all;
clc;

%% Setup
if nargin < 1
    [Type, numTestProb, params_IF, params_ID, params_CF] = setup_crossover_IF_ID_CF();
elseif nargin >=1
    [Type, numTestProb, params_IF, params_ID, params_CF] = setup_crossover_IF_ID_CF(Type);
end

% For random test only
seed = 1;                   % Seed for random number generator

% For netlib test only
nameOfProbSet = 'testNetlib.txt';

fileName = ['crossover_to_simplex_test_IF_ID_CF_' Type];

% Options for the plots
options_plot = [];
% options_evalPerf.solverNames = {'With perturbations' 'Without perturbations'};
options_plot.solverNames = {'Identification Function' 'Indicators' 'Cutoff'};
options_plot.fileName = fileName;
options_plot.logplot = 1;
options_plot.Quiet = 1;
options_plot.isCaptions = 0;

logFileName = [fileName '_log.txt'];

if exist(fullfile(cd, logFileName),'file')
    delete(logFileName);
end

diary(logFileName);
%% Run the test
fprintf('============================== Crossover Tests ==============================\n');

printHeader;

%% Initialise
switch Type
    case 'netlib'
        % read in the name of all test priblems and stoe them in a cell
        prob2test = readProbSet(nameOfProbSet);
        
        % get the number of test problems
        numTestProb = length(prob2test);
        
        %fprintf('Netlib: In total %d problems detected.\n',numTestProb);
        
    case {'random', 'random_degen'}
        rng('default');
        rng(seed);
        prob2test =  strtrim( cellstr( num2str((1:numTestProb)', 'random_%d') ) );
    otherwise
        return;
end

i = 1;

splxIter_IF = zeros(numTestProb,1);
splxIter_ID = splxIter_IF;
splxIter_CF = splxIter_IF;

mu_IF       = splxIter_IF;
mu_ID       = splxIter_IF;
mu_CF       = splxIter_IF;

ipm_iter    = splxIter_IF;
%basis_diff  = splxIter_IF;

%% Main loop
while i <= numTestProb
    switch Type
        case 'netlib'
            % load test problems
            load(prob2test{i});
            [A,b,c,FEASIBLE]=myPreprocess(A,b,c,lbounds,ubounds,BIG);
            
        case 'random'
            [A,b,c] = generateRandomProb('m_min',10,'m_max',200,...
                'n_min',20,'n_max',500);
            
        case 'random_degen'
            [A, b, c] = generateDegenProb('m_min',10,'m_max',200,...
                'n_min',20,'n_max',500);
    end
    
    %% Solve the problem using ipm with IF, ID, CF
    p_IF = pipm(A,b,c,params_IF); p_IF.solve;
    p_ID = pipm(A,b,c,params_ID); p_ID.solve;
    p_CF = pipm(A,b,c,params_CF); p_CF.solve;
    
    Prob = prob2test{i};
    
    
    %% Collect data
    splxIter_IF(i) = p_IF.getSplxIter; splxIter_ID(i) = p_ID.getSplxIter; splxIter_CF(i) = p_CF.getSplxIter;
    mu_IF(i)       = p_IF.getMu;       mu_ID(i)       = p_ID.getMu;       mu_CF(i)       = p_CF.getMu; 
    
    ipm_iter(i)     = p_IF.getIPMIterCount;
    
    printContent(Prob, p_IF, p_ID, p_CF);
    
    %% Increment counter
    i = i+1;
end
clearvars A b c lbounds ubounds NAME i Prob p_* BIG FEASIBLE ifree;
save([ 'crossover_to_simplex_test_IF_' Type '.mat']);

%% Calculate the average
fprintf('---------------------------------------------------------------------\n');
tmp_splxIter_IF = splxIter_IF; tmp_splxIter_ID = splxIter_ID; tmp_splxIter_CF = splxIter_CF;
tmp_splxIter_IF(isnan(tmp_splxIter_IF)) = [];
tmp_splxIter_ID(isnan(tmp_splxIter_ID)) = [];
tmp_splxIter_CF(isnan(tmp_splxIter_CF)) = [];
% The average value of splxIter_IF and _ID are calculated after removed
% failures.
fprintf('%10s & %4s & %4s & %9.2e & %9.2e & %9.2e & %9d & %9d & %9d & %9d\n',...
    'Average:', ' ', ' ', mean(mu_IF), mean(mu_ID), mean(mu_CF),...
    round(mean(ipm_iter)), round(mean(tmp_splxIter_IF)), round(mean(tmp_splxIter_ID)), round(mean(tmp_splxIter_CF)) );

%% Plot relative performance chart
T = [splxIter_IF splxIter_ID splxIter_CF];

% Remove problems that cannot be solved by two
indx = find(sum(isnan(T),2) > 1);
T(indx,:) = [];
fprintf('\n============================ Performance Profile ============================\n');
fprintf('\n# of Probs removed: %d\n', length(indx));
fprintf('Problems removed: \n');
fprintf('%s\n',prob2test{indx})

profiles = evalPerformance(T,options_plot);
profiles.lines = {'--' '-.' ':'};
profiles.markers= {'s' 'd' 'o'};
profiles.performaceProfile;

diary off;
end

%% %%%%% %%%%%%% %%%%%%% --- Main Func End --- %%%%%%% %%%%%%% %%%%% %%


function [Type, numTestProb, params_IF, params_ID, params_CF] = setup_crossover_IF_ID_CF(Type)
% Determine which set of problems to test on.
% Choose from the following three values:
% random, netlib, random_degen
if nargin < 1
    fprintf('Pls choose the test set [1-3]: \n');
    fprintf('\t [1]. Random test (primal nondegenerate)\n');
    fprintf('\t [2]. Random test (primal-dual degenerate)\n');
    fprintf('\t [3]. Netlib test  \n');
    usrinput_type = input('Your choice here [1-3]: ');
    if usrinput_type == 1
        Type = 'random';
    elseif usrinput_type == 2
        Type = 'random_degen';
    elseif usrinput_type == 3
        Type = 'netlib';
    else
        error('testCorrectionRatios: please choose a number from the above list');
    end
end

% IF
params_IF.verbose          = 0;
params_IF.iPer             = 0;         % no perturbations
params_IF.doCrossOver      = 1;         % conduct crossover
params_IF.mu_cap           = 1e-03;     % terminate by mu_cap and tol
params_IF.tol              = 1e-07;     % to aviod ill-conditioning
params_IF.actvPredStrtgy   = 'conservidfunc';

% ID
params_ID = params_IF;
params_ID.actvPredStrtgy   = 'conservindica';

% Cutoff
params_CF = params_IF;
params_CF.actvPredStrtgy   = 'conservcutoff';

numTestProb = 100;          % Set to 10 for demo. 100 for real test.
end

function printHeader
% Header: 1       2   3      4     5     6     7   8     9
fprintf('%10s & %4s & %4s & %9s & %9s & %9s & %9s & %9s & %9s & %9s \\\\ \n',...
    'Prob', 'm', 'n', 'mu_IF', 'mu_ID', 'mu_CF', 'iter_ipm', 'splx_IF', 'splx_ID', 'splx_CF');
end

function printContent(Prob, p_IF, p_ID, p_CF)
% Iter:  1       2     3     4       5       6      7      8     9       10
fprintf('%10s & %4d & %4d & %9.2e & %9.2e & %9.2e & %9d & %9d & %9d & %9d \\\\ \n',...
    Prob, p_IF.prob.m, p_IF.prob.n, p_IF.getMu, p_ID.getMu, p_CF.getMu,...
    p_IF.getIPMIterCount, p_IF.getSplxIter, p_ID.getSplxIter, p_CF.getSplxIter);
end
