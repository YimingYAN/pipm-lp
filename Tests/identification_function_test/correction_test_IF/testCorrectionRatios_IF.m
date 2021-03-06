function testCorrectionRatios_IF
% FUNCTION testCorrectionRatios_IF
% Run correction ratios tests wtih identification function
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
% 3. Predicted original active-set from unperturbed alg \
%                           V.S.                        | Alg 6.2 - Splx
%    Actual perturbed activ-set from linprog (simplex)  /
%
% 4. Predicted original active-set from unperturbed alg \
%                           V.S.                        | Alg 6.2 - IPM
%    Actual perturbed activ-set from linprog (ipm)      /
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
%
%
% Date        :  24 March 2014
% Author      :  Yiming Yan
% Affiliation :  University of Edinburgh

close all;
clc;

%% Setup
[Type, numTestProb, params_per, params_unper] = setup_correctionRatio_IF;

% -------------------------------------------------------------------------
stopAtRangeL = 6;
stopAtRangeU = 18;

% Options for plots
Legends = { 'Alg 6.1 - Splx'  'Alg 6.1 - IPM'  'Alg 6.2 - Splx' 'Alg 6.2 - IPM' };
colors = {'r' 'b' 'k' [0 .5 0]};
lineStyles = {'-' '--' '-' '--'};
markers = {'o' '*' 's' 'd'};
fileName = ['correction_ratio_test_IF_' Type];

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
        
        %% Predict the original actv using pipm-lp (perturbed alg)
        params_per.maxIter = k;
        per = pipm(A,b,c,params_per); per.solve;
        
        %% Predict the original actv using pipm-lp iPer = 0 (unperturbed alg)
        params_unper.maxIter = k;
        unper = pipm(A, b, c, params_unper); unper.solve;
         
        %% Get the correction ratios
        
        % Predicted original actv from pipm-lp VS original actv from simplex
        [fpr_per_splx(i), mpr_per_splx(i), cr_per_splx(i)] = ...
            getCorrectionRatio(per.getActv, actualActv_splx);
        
        % Predicted original actv from pipm-lp VS original actv from ipm
        [fpr_per_ipm(i), mpr_per_ipm(i), cor_per_ipm(i)] = ...
            getCorrectionRatio(per.getActv, actualActv_ipm);
        
        % Predicted original actv from pipm-lp with iPer=0 VS original actv 
        % from simplex
        [fpr_unp_splx(i), mpr_unp_splx(i), cor_unp_splx(i)] = ...
            getCorrectionRatio(unper.getActv, actualActv_splx);
        
        % Predicted original actv from pipm-lp with iPer=0 VS original actv
        % from ipm
        [fpr_unp_ipm(i), mpr_unp_ipm(i), cor_unp_ipm(i)] = ...
            getCorrectionRatio(unper.getActv, actualActv_ipm);
        
        avgRes_per(i)   = per.getIPMResidual;
        avgRes_unper(i) = unper.getIPMResidual;
        
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
    
    %save
    
    counter = counter+1;
end

fprintf('\n\tTotal number of probs skipped: %d\n', skipped);

clearvars A b c lb ub m n i k xsol exitflag...
    Prob per unper cr* fpr* mpr* tp*...
    avgRes1 avgRes2 actualActv counter;

%% Output the result
fprintf('\n3. Output the result...\n');
range = stopAtRangeL : stopAtRangeU;

save( ['correction_ratio_test_IF_' Type '.mat'] );

plotCorrectionRatios_4lines(falsePrediction, missedPrediction,...
    correctionR, avgResidual, range, Legends, fileName,...
    colors, lineStyles, markers);

fprintf('DONE.\n');
fprintf('Pls check the file %s for the plot.\n', [fileName '.pdf']);
end

%% ----------------- Main Func End ----------------- %%

%% Print iterative info
function printHeader
fprintf('\n%4s | %11s | %7s %7s %7s %9s | %7s %7s %7s %9s | %7s %7s %7s %9s | %7s %7s %7s %9s\n',...
    'Iter', '[m,n]  ',...
    'F_Per_S', 'M_Per_S', 'C_Per_S', 'R_Per_S',...
    'F_Per_I', 'M_Per_I', 'C_Per_I', 'R_Per_I',...
    'F_Unp_S', 'M_Unp_S', 'C_Unp_S', 'R_Unp_S',...
    'F_Unp_I', 'M_Unp_I', 'C_Unp_I', 'R_Unp_I');
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

function plotCorrectionRatios_4lines(falsePrediction, missedPrediction,...
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

