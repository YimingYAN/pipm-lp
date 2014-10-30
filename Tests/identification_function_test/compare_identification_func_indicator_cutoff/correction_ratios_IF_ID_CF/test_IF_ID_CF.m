function test_IF_ID_CF(Type)
% test_IF_ID_CF This function compares the accuracy of predicted active-sets
% from primal-dual pathfollowing with identification function, indicators,
% and simple cutoff.
%
% March 25, 2014
% Yiming Yan

%% Determine which set of problems to test on.
% Choose from the following two values:
% random, random_degen, netlib
if nargin < 1
    fprintf('1. Pls choose the test set [1-2]: \n');
    fprintf('\t [1]. Random test (primal nondegenerate)\n');
    fprintf('\t [2]. Random test (primal-dual degenerate)\n');
    fprintf('\t [3]. 6 netlib problems\n')
    usrinput_type = input('Your choice here [1-2]: ');
    
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


%% Set Parameters for ipm
nameOfProbSet = 'netlib_IF_ID_CF.txt';
% nameOfProbSet = 'netlib_crossover_IF_ID_CF.txt';

stopAtRangeL = 6;
stopAtRangeU = 18;

% IF
params_IF.verbose          = 0;
params_IF.iPer             = 0;         % no perturbations
params_IF.doCrossOver      = 0;
params_IF.mu_cap           = 1e-07;     % terminate by mu_cap and tol
params_IF.tol              = 1e-07;     % to aviod ill-conditioning
params_IF.actvPredStrtgy   = 'conservidfunc';

% ID
params_ID = params_IF;
params_ID.actvPredStrtgy   = 'conservindica';

% Cutoff
params_CF = params_IF;
params_CF.actvPredStrtgy   = 'conservcutoff';

% Options for plots
Legends = { 'Identification function'  'Indicators'  'Cut-off' };
colors = {'r' 'b' 'k'};
lineStyles = {'-' '--' '-.'};
markers = {'o' '*' 's'};
fileName = ['correction_ratio_test_IF_ID_CF_' Type];

%% Deterine test problems
switch Type
    case {'random', 'random_degen'}
        numTestProb  = 100; % set to 10 for demo. 100 for real test.
    case 'netlib'
        % read in the name of all test priblems and stoe them in a cell
        prob2test = readProbSet(nameOfProbSet);
        
        % get the number of test problems
        numTestProb = length(prob2test);
end


%% Run the test
logFileName = [fileName '_log.txt'];
if exist(logFileName, 'file')
    delete(logFileName);
end
diary(logFileName);
fprintf('\n2. Start the %s test...\n', Type);
fprintf('\n================================= Correction Ratio Tests =================================\n');
printHeader;

counter = 1;      % Counter for outter loop

falsePrediction  = zeros(stopAtRangeU - stopAtRangeL, 3);
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
    
    fpr_IF = zeros(numTestProb,1); mpr_IF = fpr_IF; cr_IF  = fpr_IF;
    
    fpr_ID = fpr_IF; mpr_ID = fpr_IF; cr_ID = fpr_IF;
    
    fpr_CF  = fpr_IF; mpr_CF  = fpr_IF; cr_CF  = fpr_IF;

    
    avgRes_IF = fpr_IF; avgRes_ID = fpr_IF; avgRes_CF = fpr_IF;
    
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
            case 'netlib'
                load(prob2test{i});
                 [A,b,c,FEASIBLE]=myPreprocess(A,b,c,lbounds,ubounds,BIG);
            otherwise
                return;
        end
        
        
        %% Solve for 'actual active-set' of the original problem
        [m, n] = size(A);
                
        %  Get the actual original actv from linprog (ipm)
        [ actualActv, exitflag  ] = solveLinprog(A, b, c,'ipm' );
        
        % Skip the test problem if linprog does not converge
        if exitflag ~= 1
            skipped = skipped + 1;
            if strcmpi(Type, 'netlib')
                i=i+1;
            end
            continue;
        end
        
        %% Predict the actv - IF
        params_IF.maxIter = k;
        p_IF = pipm(A,b,c,params_IF); p_IF.solve;
        
        %% Predict the actv - ID
        params_ID.maxIter = k;
        p_ID = pipm(A, b, c, params_ID); p_ID.solve;
        
        %% Predict the actv - CF
        params_CF.maxIter = k;
        p_CF = pipm(A,b,c,params_CF); p_CF.solve;
        
        %% Get the correction ratios
        % Actv - IF
        [fpr_IF(i), mpr_IF(i), cr_IF(i)] = ...
            getCorrectionRatio(p_IF.getActv, actualActv);
        
        % Actv - ID
        [fpr_ID(i), mpr_ID(i), cr_ID(i)] = ...
            getCorrectionRatio(p_ID.getActv, actualActv);
        
        % Actv - CF
        [fpr_CF(i), mpr_CF(i), cr_CF(i)] = ...
            getCorrectionRatio(p_CF.getActv, actualActv);
     
        avgRes_IF(i)   = p_IF.getIPMResidual;
        avgRes_ID(i)   = p_ID.getIPMResidual;
        avgRes_CF(i)   = p_CF.getIPMResidual;
        
        Avgm(i) = m; Avgn(i) = n;
        
        %% Increment the counter
        i = i+1;
    end
    
    %% Get the averages
    falsePrediction(counter,:)  = mean( [ fpr_IF  fpr_ID fpr_CF] );
    missedPrediction(counter,:) = mean( [ mpr_IF  mpr_ID mpr_CF] );
    correctionR(counter,:)      = mean( [ cr_IF   cr_ID  cr_CF ] );
    
    avgResidual(counter,:)      = mean([avgRes_IF avgRes_ID avgRes_CF]);
    
    Avgm = round(mean(Avgm));   Avgn = round(mean(Avgn));
    
    printContent(k, counter, Avgm, Avgn, falsePrediction,...
        missedPrediction, correctionR, avgResidual);
    
    %save
    
    counter = counter+1;
end

fprintf('\n\tTotal number of probs skipped: %d\n', skipped);

clearvars A b c lb ub m n i k xsol exitflag...
    Prob per unper actualActv counter;

%% Output the result
fprintf('\n3. Output the result...\n');
range = stopAtRangeL : stopAtRangeU;

save( [fileName '.mat'] );

plotCorrectionRatios_3lines(falsePrediction, missedPrediction,...
    correctionR, avgResidual, range, Legends, fileName,...
    colors, lineStyles, markers);

fprintf('DONE.\n');
fprintf('Pls check the file %s for the plot.\n', [fileName '.pdf']);
diary off;
end

%% Print iterative info
function printHeader
fprintf('\n%4s | %11s | %7s %7s %7s %9s | %7s %7s %7s %9s | %7s %7s %7s %9s\n',...
    'Iter', '[m,n]  ',...
    'F_IF', 'M_IF', 'C_IF', 'R_IF',...
    'F_ID', 'M_ID', 'C_ID', 'R_ID',...
    'F__CF', 'M__CF', 'C__CF', 'R__CF');
end

function printContent(k, counter, Avgm, Avgn, falsePrediction,...
    missedPrediction, correctionR, avgResidual)
fprintf('%4d | [%4d %4d] | %7.2f %7.2f %7.2f %9.2e | %7.2f %7.2f %7.2f %9.2e | %7.2f %7.2f %7.2f %9.2e\n',...
    k,Avgm,Avgn,...
    falsePrediction(counter, 1), missedPrediction(counter, 1), ...
    correctionR(counter, 1),      avgResidual(counter, 1), ...
    falsePrediction(counter, 2), missedPrediction(counter, 2), ...
    correctionR(counter, 2),      avgResidual(counter, 2),...
    falsePrediction(counter, 3), missedPrediction(counter, 3), ...
    correctionR(counter, 3),      avgResidual(counter, 3));

end

function plotCorrectionRatios_3lines(falsePrediction, missedPrediction,...
    correctionR, avgResidual, range, Legends, fileName,...
    colors, lineStyles, markers)
clf;

titles = {'False Prediction Ratio',...
    'Missed Prediction Ratio',...
    'Correction Ratio',...
    'Average Relative Residual (log10)'};

h = zeros(4,1);
for i = 1: 4
    h(i) = subplot(2,2,i);
    
    switch i
        case 1, T = falsePrediction;
        case 2, T = missedPrediction;
        case 3, T = correctionR;
        case 4, T = log10(avgResidual);
    end
    
    hold on
    for j = 1 : size(T,2)
        plot(range,T(:,j),...
            [lineStyles{j} markers{j}],...
            'Color', colors{j},...
            'MarkerSize', 8);%,...
            %'MarkerFaceColor', colors{j})
    end
    hold off
    
    title(titles{i});
    
end

% Set axes properties
set(h,'XTick',range);
set(h,'XLim',[range(1) range(end)+0.1]);
set(h(1:3),'YLim',[-0.1 1.1]);
set(h(4), 'YLim',[floor(min(min(log10(avgResidual)))) ceil(max(max(log10(avgResidual))))]);
set(h,'XGrid','on','YGrid','on');
set(h(4),'YTick', floor(min(min(log10(avgResidual)))):1:ceil(max(max(log10(avgResidual)))));
hleg = legend(Legends,'Orientation','horizontal');
p =get(hleg,'Position');
p(1) = 0.5-0.48*p(3); p(2)= 0.02;
set(hleg,'Position',p);
print('-depsc',fileName);
try
    [result,msg] = eps2pdf([fileName '.eps']);
    if result
        fprintf('epstopdf failed. No pdf file generated.\n');
    end
catch err_plotcr
    fprintf('epstopdf failed. No pdf file generated.\n');
end

end