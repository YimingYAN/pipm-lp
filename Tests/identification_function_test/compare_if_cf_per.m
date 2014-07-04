% This script is used to compare the performance of perturbed alg with
% identificaiton fucntion and perturbed alg with cutoff

clear;
clc;
type = {'random', 'random_degen'};

T = [];
% T = [ perAlg_CF_vs_splx unperAlg_CF_vs_ipm perAlg_IF_vs_ipm
% unperAlg_IF_vs_ipm ];

for i = 1:2
    
    %% load cr test CF random
    load(['../correction_ratio_test/correction_ratio_test_' type{i} '.mat']);
    
    T(:,1:2)  = [ correctionR(:,1) correctionR(:,4) ];
    TR(:,1:2) = log10( [ avgResidual(:,1) avgResidual(:,4) ] );
    
    %% load cr test IF random
    load(['correction_test_IF/correction_ratio_test_' type{i} '.mat']);
    
    T(:,3:4) = [ correctionR(:,2) correctionR(:,4) ];
    TR(:,3:4) = log10( [ avgResidual(:,2) avgResidual(:,4) ] );
    
    %% plot
    figure;
    titles = {'Correction Ratio',...
    'Average Relative Residual (log10)'};
    fileName = ['compare_if_cf_per_' type{i}];
    Legends = { 'Alg 6.1 (CF) - Splx'  'Alg 6.2 (CF) - IPM'  'Alg 6.1 (IF) - IPM' 'Alg 6.2 (IF) - IPM' };
    colors = {'r' 'b' 'k' [0 .5 0]};
    lineStyles = {'-' '--' '-' '--'};
    markers = {'o' '*' 's' 'd'};
    
    h(1) = subplot(1,2,1);
    hold on
    for j = 1 : size(T,2)
        plot(range,T(:,j),...
            [lineStyles{j} markers{j}],...
            'Color', colors{j},...
            'MarkerSize', 8);
    end
    hold off
    
    title(titles{1});
    
    
    h(2) = subplot(1,2,2);
    hold on
    for j = 1 : size(T,2)
        plot(range,TR(:,j),...
            [lineStyles{j} markers{j}],...
            'Color', colors{j},...
            'MarkerSize', 8);
    end
    
    hold off
    title(titles{2});
    
    % Set figure properties
    % Set axes properties
    axis(h, 'square');
    set(h,'XTick',range);
    set(h,'XLim',[range(1) range(end)]);
    set(h(2),'YTick', floor(min(min(log10(avgResidual)))):1:ceil(max(max(log10(avgResidual)))));
    set(h(1),'YLim',[0 1]);
    set(h,'XGrid','on','YGrid','on');
    
    hleg = legend(Legends,'Orientation','horizontal');
    set(hleg,'FontSize',9);
    p =get(hleg,'Position');
    
    p(1) = 0.5-0.48*p(3); p(2)= 0.2; 
    set(hleg,'Position',p);
    

    print('-depsc',fileName);
end

