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