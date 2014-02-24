% This script combines all cr plots
function combine_cr_plots
range = 6:18; 
L = range(end) - range(1) + 1;

falsePrediction_4lines = [];
missedPrediction_4lines = [];
correctionR_4lines = [];
avgResidual_4lines = [];

% random
% per vs splx; unper vs ipm
load correction_ratio_test_random_per_splx_unper_ipm.mat

falsePrediction_4lines = [falsePrediction_4lines falsePrediction];
missedPrediction_4lines = [missedPrediction_4lines missedPrediction];
correctionR_4lines = [correctionR_4lines correctionR];
avgResidual_4lines = [avgResidual_4lines avgResidual];

% per vs ipm;  unper vs splx
load correction_ratio_test_random_per_ipm_unper_splx.mat

falsePrediction_4lines = [falsePrediction_4lines falsePrediction];
missedPrediction_4lines = [missedPrediction_4lines missedPrediction];
correctionR_4lines = [correctionR_4lines correctionR];
avgResidual_4lines = [avgResidual_4lines avgResidual];

falsePrediction_4lines(L+1:end,:) = []; 
missedPrediction_4lines(L+1:end,:) = []; 
correctionR_4lines(L+1:end,:) = []; 
avgResidual_4lines(L+1:end,:) = [];

fileName = [fileName '_4lines'];
Legends = {'Alg. 6.1 - Splx' 'Alg. 6.2 - IPM' 'Alg. 6.1 - IPM' 'Alg. 6.2 - Splx'};

plotCorrectionRatios_4lines(falsePrediction_4lines, missedPrediction_4lines,...
    correctionR_4lines, avgResidual_4lines, range, Legends, fileName);

pause


% random degenerate
falsePrediction_4lines = [];
missedPrediction_4lines = [];
correctionR_4lines = [];
avgResidual_4lines = [];

% per vs splx; unper vs ipm
load correction_ratio_test_random_degen_per_splx_unper_ipm_20.mat

falsePrediction_4lines = [falsePrediction_4lines falsePrediction];
missedPrediction_4lines = [missedPrediction_4lines missedPrediction];
correctionR_4lines = [correctionR_4lines correctionR];
avgResidual_4lines = [avgResidual_4lines avgResidual];

% per vs ipm_unper vs splx
load correction_ratio_test_random_degen_per_ipm_unper_splx_20.mat

falsePrediction_4lines = [falsePrediction_4lines falsePrediction];
missedPrediction_4lines = [missedPrediction_4lines missedPrediction];
correctionR_4lines = [correctionR_4lines correctionR];
avgResidual_4lines = [avgResidual_4lines avgResidual];

falsePrediction_4lines(L+1:end,:) = []; 
missedPrediction_4lines(L+1:end,:) = []; 
correctionR_4lines(L+1:end,:) = []; 
avgResidual_4lines(L+1:end,:) = [];

fileName = [fileName '_4lines'];
Legends = {'Alg. 6.1 - Splx' 'Alg. 6.2 - IPM' 'Alg. 6.1 - IPM' 'Alg. 6.2 - Splx'};

plotCorrectionRatios_4lines(falsePrediction_4lines, missedPrediction_4lines,...
    correctionR_4lines, avgResidual_4lines, range, Legends, fileName);
end

%% This function is used to plot the correction ratios
function plotCorrectionRatios_4lines(falsePrediction, missedPrediction,...
    correctionR, avgResidual, range, Legends, fileName)
clf;
lineStyle = {'-ro' '--b*' '-k+' '--mv'};

h(1) = subplot(2,2,1);
plot(range,falsePrediction(:,1),lineStyle{1},...
    range,falsePrediction(:,2),lineStyle{2},...
    range,falsePrediction(:,3),lineStyle{3},...
    range,falsePrediction(:,4),lineStyle{4});
title('False Prediction Ratio');

h(2) = subplot(2,2,2);
plot(range,missedPrediction(:,1),lineStyle{1},...
    range,missedPrediction(:,2),lineStyle{2},...
    range,missedPrediction(:,3),lineStyle{3},...
    range,missedPrediction(:,4),lineStyle{4});
title('Missed Prediction Ratio');

h(3) = subplot(2,2,3);
plot(range,correctionR(:,1),lineStyle{1},...
    range,correctionR(:,2),lineStyle{2},...
    range,correctionR(:,3),lineStyle{3},...
    range,correctionR(:,4),lineStyle{4});
title('Correction Ratio');

h(4) = subplot(2,2,4);
plot(range,log10(avgResidual(:,1)),lineStyle{1},...
    range,log10(avgResidual(:,2)),lineStyle{2},...
    range,log10(avgResidual(:,3)),lineStyle{3},...
    range,log10(avgResidual(:,4)),lineStyle{4})
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
p(1) = 0.5-0.48*p(3); p(2)= 0.02;
set(hleg,'Position',p);
print('-depsc',fileName);
[result,msg] = eps2pdf([fileName '.eps']);

end

