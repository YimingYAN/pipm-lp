function testCR_IF_ID_CF_Netlib_M
% testCR_IF_ID_CF_Netlib_M This function plots the prediction ratios for 
% netlib problems. 
% 1. We solve the problem to optimality and record the total number of ipm iterations, M
% 2. Run the algorithm again and compare the prediction rations at the last
% several itertaoins for each problem.
%
% 21 June 2014
% Yiming Yan

%% Determine which set of problems to test on.
Type = 'netlib';


%% Set Parameters for ipm
% nameOfProbSet = 'netlib_IF_ID_CF.txt';
nameOfProbSet = 'netlib_crossover_IF_ID_CF.txt';

steps = 6;

% Per
params.verbose          = 0;
params.iPer             = 0;         % no perturbations
params.doCrossOver      = 0;
params.mu_cap           = 1e-06;     % terminate by mu_cap and tol
params.tol              = 1e-06;     % to aviod ill-conditioning

% IF
params_IF = params;
params_IF.actvPredStrtgy   = 'conservidfunc';

% ID
params_ID = params;
params_ID.actvPredStrtgy   = 'conservindica';

% Cutoff
params_CF = params;
params_CF.actvPredStrtgy   = 'conservcutoff';


% Options for plots
Legends = { 'Identification function'  'Indicators'  'Cut-off' };
colors = {'r' 'b' 'k'};
lineStyles = {'-' '--' '-.'};
markers = {'o' '*' 's'};
fileName = ['correction_ratio_test_IF_ID_CF_' Type '_Netlib_M'];

%% Deterine test problems

% read in the name of all test priblems and stoe them in a cell
prob2test = readProbSet(nameOfProbSet);

% get the number of test problems
numTestProb = length(prob2test);

%% Run the test
logFileName = [fileName '_log.txt'];
if exist(logFileName, 'file')
    delete(logFileName);
end
diary(logFileName);
fprintf('\n2. Start the %s test...\n', Type);
fprintf('\n================================= Correction Ratio Tests =================================\n');

i = 1;

fprintf('Progress: \n');
while i <= numTestProb
    fprintf('| --> %d',i);
    
    load(prob2test{i});
    [A,b,c,FEASIBLE]=myPreprocess(A,b,c,lbounds,ubounds,BIG);
    
    %% Solve for 'actual active-set' of the original problem
    [m, n] = size(A);
    
    %  Get the actual original actv from linprog (ipm)
    [ actualActv, exitflag  ] = solveLinprog(A, b, c,'ipm' );
    
    % Skip the test problem if linprog does not converge
    if exitflag ~= 1
        skipped = skipped + 1;
        i=i+1;
        continue;
    end
    
    % Determine range
    params.maxIter = 100;
    p = pipm(A,b,c,params); p.solve;
    stopAtRangeU = p.getIPMIterCount;
    stopAtRangeL = stopAtRangeU - steps + 1;
    
    counter = 1;      % Counter for outter loop
    skipped = 0; 
    for k = stopAtRangeL:1:stopAtRangeU
        
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
        [fpr_IF(counter, i), mpr_IF(counter, i), cr_IF(counter, i)] = ...
            getCorrectionRatio(p_IF.getActv, actualActv);
        
        % Actv - ID
        [fpr_ID(counter, i), mpr_ID(counter, i), cr_ID(counter, i)] = ...
            getCorrectionRatio(p_ID.getActv, actualActv);
        
        % Actv - CF
        [fpr_CF(counter, i), mpr_CF(counter, i), cr_CF(counter, i)] = ...
            getCorrectionRatio(p_CF.getActv, actualActv);
        
        res_IF(counter, i) = p_IF.getIPMResidual
        res_ID(counter, i) = p_ID.getIPMResidual
        res_CF(counter, i) = p_CF.getIPMResidual
        
        counter = counter+1;
    end
    fprintf(' |\n');
    %% Increment the counter
    i = i+1;
end

%% Get the matrix
falsePrediction  = [ mean(fpr_IF,2)  mean(fpr_ID,2) mean(fpr_CF,2)];
missedPrediction = [ mean(mpr_IF,2)  mean(mpr_ID,2) mean(mpr_CF,2)];
correctionR      = [ mean(cr_IF,2)   mean(cr_ID,2)  mean(cr_CF,2) ];
Residual         = [ mean(res_IF,2)  mean(res_ID,2) mean(res_CF,2)];

fprintf('\n\tTotal number of probs skipped: %d\n', skipped);

clearvars A b c lb ub m n i k xsol exitflag...
    Prob per unper actualActv counter;

%% Output the result
fprintf('\n3. Output the result...\n');
range = 1 : steps;

save( [fileName '.mat'] );

plotCorrectionRatios_3lines(falsePrediction, missedPrediction,...
    correctionR, Residual, range, Legends, fileName,...
    colors, lineStyles, markers);

fprintf('DONE.\n');
fprintf('Pls check the file %s for the plot.\n', [fileName '.pdf']);
diary off;
end

%% Print iterative info
function plotCorrectionRatios_3lines(falsePrediction, missedPrediction,...
    correctionR, Residual, range, Legends, fileName,...
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
        case 4, T = log10(Residual);
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
for i = 1:length(range)
    if i ~= 1
        tickLabel{length(range) - i + 1} = ['M-' num2str(range(i)) - 1] ;
    else
        tickLabel{length(range) - i + 1} = ['M'] ;
    end
end
set(h,'XTick',range);
set(h,'XTickLabel', tickLabel);
set(h,'XLim',[range(1) range(end)+0.1]);
set(h(1:3),'YLim',[-0.1 1.1]);
set(h(4), 'YLim',[floor(min(min(log10(Residual)))) ceil(max(max(log10(Residual))))]);
set(h,'XGrid','on','YGrid','on');
set(h(4),'YTick', floor(min(min(log10(Residual)))):1:ceil(max(max(log10(Residual)))));
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
