function testCorrectionRatios_Netlib
% testCorrectionRatios_Netlib This function plots the prediction ratios for
% netlib problems.
% 1. We solve the problem to optimality and record the total number of ipm iterations, M
% 2. Run the algorithm again and compare the prediction rations at the last
% several itertaoins for each problem.
%
% 21 June 2014
% Yiming Yan

%% Determine which set of problems to test on.
Type = 'netlib';
warning off

%% Set Parameters for ipm
% nameOfProbSet = 'netlib_IF_ID_CF.txt';
nameOfProbSet = 'netlib_test_probs.txt';

steps = 10;

actvPredStrtgy = 'conservCutoff';

% Without perturbations
parameters_unper.verbose        = 0;
parameters_unper.iPer           = 0;
parameters_unper.actvPredStrtgy = actvPredStrtgy;
parameters_unper.doCrossOver    = 0;
parameters_unper.mu_cap         = 1e-08;     
parameters_unper.tol            = 1e-08;     

% With perturbations
parameters_per.verbose          = 0;
parameters_per.iPer             = 1e-02;
parameters_per.actvPredStrtgy   = actvPredStrtgy;
parameters_per.doCrossOver      = 0;
parameters_per.mu_cap           = 1e-08;     
parameters_per.tol              = 1e-08;     

% Ref
parameters_ref = parameters_unper;
parameters_ref.maxIter = 100;

% Options for plots
Legends = { 'Alg 6.1 - Splx'  'Alg 6.1 - IPM'  'Alg 6.2 - Splx' 'Alg 6.2 - IPM' };
colors = {'r' 'b' 'k' [0 .5 0]};
lineStyles = {'-' '--' '-' '--'};
markers = {'o' '*' 's' 'd'};
fileName = 'testCorrectionRatios_Netlib';

%% Determine test problems

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
fprintf('\n======================= Correction Ratio Tests (Netlib) =======================\n');

i = 1;

fpr_per_splx = NaN*ones(steps,numTestProb);
mpr_per_splx = fpr_per_splx;
cor_per_splx = fpr_per_splx;
fpr_per_ipm  = fpr_per_splx;
mpr_per_ipm  = fpr_per_splx;
cor_per_ipm  = fpr_per_splx;
fpr_unp_splx = fpr_per_splx;
mpr_unp_splx = fpr_per_splx;
cor_unp_splx = fpr_per_splx;
fpr_unp_ipm  = fpr_per_splx;
mpr_unp_ipm  = fpr_per_splx;
cor_unp_ipm  = fpr_per_splx;
res_per = fpr_per_splx;
res_unp = fpr_per_splx;

fprintf('Progress: \n');
while i <= numTestProb
    fprintf('| --> %2d',i);
    
    load(prob2test{i});
    [A,b,c,~]=myPreprocess(A,b,c,lbounds,ubounds,BIG);
    
    %% Solve for 'actual active-set' of the original problem    
    %  Get the actual original actv from linprog (simplex)
    [ actualActv_splx, exitflag_splx ] = solveLinprog(A, b, c,'splx');
    
    %  Get the actual original actv from linprog (ipm)
    [ actualActv_ipm,  exitflag_ipm  ] = solveLinprog(A, b, c,'ipm');
    
    % Skip the test problem if linprog does not converge
    if exitflag_splx ~= 1 || exitflag_ipm ~= 1
        skipped = skipped + 1;
        i = i + 1;
        fprintf(' (-) skipped due to linprog \n');
        continue;
    end
    
    % Determine range
    ref = pipm(A,b,c,parameters_ref); ref.solve;
    stopAtRangeU = ref.getIPMIterCount;
    stopAtRangeL = stopAtRangeU - steps + 1;
    
    counter = 1;      % Counter for outter loop
    skipped = 0;
    
    for k = stopAtRangeL:1:stopAtRangeU
        
        %% Predict the actv - Per
        parameters_per.maxIter = k;
        per = pipm(A,b,c,parameters_per); per.solve;
        
        %% Predict the actv - Unper
        parameters_unper.maxIter = k;
        unper = pipm(A, b, c, parameters_unper); unper.solve;
        
        %% Get the correction ratios
        % Predicted original actv from pipm-lp VS original actv from simplex
        [fpr_per_splx(counter, i), mpr_per_splx(counter, i), cor_per_splx(counter, i)] = ...
            getCorrectionRatio(per.getActv, actualActv_splx);
        
        % Predicted original actv from pipm-lp VS original actv from ipm
        [fpr_per_ipm(counter, i), mpr_per_ipm(counter, i), cor_per_ipm(counter, i)] = ...
            getCorrectionRatio(per.getActv, actualActv_ipm);
        
        % Predicted original actv from pipm-lp with iPer=0 VS original actv
        % from simplex
        [fpr_unp_splx(counter, i), mpr_unp_splx(counter, i), cor_unp_splx(counter, i)] = ...
            getCorrectionRatio(unper.getActv, actualActv_splx);
        
        % Predicted original actv from pipm-lp with iPer=0 VS original actv
        % from ipm
        [fpr_unp_ipm(counter, i), mpr_unp_ipm(counter, i), cor_unp_ipm(counter, i)] = ...
            getCorrectionRatio(unper.getActv, actualActv_ipm);
        
        res_per(counter, i) = per.getIPMResidual;
        res_unp(counter, i) = unper.getIPMResidual;
        
        counter = counter+1;
    end
    
    fprintf(' |\n');
    %% Increment the counter
    i = i+1;
end

fprintf('\n');

%% Remove problems that cannot solved
idx_per = sum(isnan(res_per),1) == 0;
idx_unp = sum(isnan(res_unp),1) == 0;
idx = idx_per + idx_unp > 1;

%% Get the matrix
falsePrediction  = [ mean(fpr_per_splx(:,idx),2) mean(fpr_per_ipm(:,idx),2) mean(fpr_unp_splx(:,idx),2) mean(fpr_unp_ipm(:,idx),2)]
missedPrediction = [ mean(mpr_per_splx(:,idx),2) mean(mpr_per_ipm(:,idx),2) mean(mpr_unp_splx(:,idx),2) mean(mpr_unp_ipm(:,idx),2)]
correctionR      = [ mean(cor_per_splx(:,idx),2) mean(cor_per_ipm(:,idx),2) mean(cor_unp_splx(:,idx),2) mean(cor_unp_ipm(:,idx),2)]

avgResidual      = [ mean(res_per(:,idx),2) mean(res_per(:,idx),2) mean(res_unp(:,idx),2) mean(res_unp(:,idx),2)]

fprintf('\n\tTotal number of probs skipped: %d\n', skipped);

clearvars A b c lb ub m n i k xsol exitflag*...
    Prob per unper actualActv* counter...
    NAME ans ;

%% Output the result
fprintf('\n3. Output the result...\n');
range = 1 : steps;

save( [fileName '.mat'] );

plotCorrectionRatios_4lines(falsePrediction, missedPrediction,...
    correctionR, avgResidual, range, Legends, fileName,...
    colors, lineStyles, markers, 1);

fprintf('DONE.\n');
fprintf('Pls check the file %s for the plot.\n', [fileName '.pdf']);
diary off;
end