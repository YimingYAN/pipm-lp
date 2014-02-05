% This script is used to check the initial value of mu.
clear all;
close all;
clc;

%% Setup
% Determine which set of problems to test on.
% Choose from the following three values:
% random, netlib, random_degen
fprintf('1. Pls choose the test set [1-3]: \n');
fprintf('\t 1. Random test (primal nondegenerate)\n');
fprintf('\t 2. Random test (primal-dual degenerate)\n');
fprintf('\t 3. Netlib test  \n');
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

% -------------------------------------------------------------------------

% For random test only
numTestProb = 100;          % Set to 10 for demo. 100 for real test.
seed = 1;                   % Seed for random number generator


% For netlib test only
nameOfProbSet = 'testNetlib.txt';

% Set initial perturbations
iPer = 1e-02;

%% Run the test
diary(['checkInitialMu_' Type '.log'])
fprintf('=============== Compare mu0 ===============\n');

% Initialize
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

mup_0 = 0; mu_0 = 0;
mus = zeros(numTestProb,2);
i = 1;
fprintf('%10s %4s %9s %9s %9s\n', 'ID.', 'N', 'MUP_0', 'MU_0', 'REL_DIFF');
while i<=numTestProb
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
    
    % Get (x0,y0,s0)
    prob = Prob(full(A),full(b),full(c));
    iter = Iterate(prob);
    iter.initialPoint(prob);
    
    % Get mup_0 and mu_0
    mup_0 = mean((prob.x+iPer*ones(prob.n,1) ).*(prob.s+iPer*ones(prob.n,1)));
    mu_0  = mean( prob.x.*prob.s);
    
    mus(i,:) = [ mup_0 mu_0 ];
    
    fprintf('%10s %4d %9.2e %9.2e %9.2e\n',...
        prob2test{i}, prob.n, mup_0, mu_0, (mup_0-mu_0)/mu_0);
    
    % Increment counter
    i = i+1;
end
fprintf('----------------------------------------------\n')
fprintf('%10s %4s %9.2e %9.2e %9.2e\n',...
    'MEAN', '--', mean(mus(:,1)), mean(mus(:,2)),...
    mean( ( abs(mus(:,1) - mus(:,2)) ) ./ mus(:,2)));

diary off;
