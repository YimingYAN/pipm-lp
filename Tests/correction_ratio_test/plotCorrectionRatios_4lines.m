%% Function used to calculate the correction ratios - 4 lines
function plotCorrectionRatios_4lines(falsePrediction, missedPrediction,...
    correctionR, avgResidual, range, Legends, fileName, colors,...
    lineStyles, markers, isNetlib)
% PLOTCORRECTIONRATIOS This function plots the prediction ratios, including 
% false prediction, missed prediction and correction ratios.
%
% Date  : 22 Feb 2014
% Author: Yiming Yan

clf;

if nargin < 11
    isNetlib = 0;
end
titles = {'False Prediction Ratio',...
    'Missed Prediction Ratio',...
    'Correction Ratio',...
    'Average Relative Residual (log10)'};

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontSize',9);
set(0,'DefaultTextFontSize',9);
set(0,'defaulttextinterpreter','tex');

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
    title(titles{i}, 'FontName', 'Arial', 'FontWeight', 'Normal');
end

% Set axes properties
set(h,'XTick',range);
if isNetlib
    for i = 1:length(range)
        if i ~= 1
            tickLabel{length(range) - i + 1} = ['M-' num2str(range(i) - 1)];
        else
            tickLabel{length(range) - i + 1} = ['M'];
        end
    end
    set(h,'XTickLabel', tickLabel);
end
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