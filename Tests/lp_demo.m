% lp_demo - a demo for comparing the perturbed and unperturbed algorithm.
%     
% iteration information for both algorithms are obtained and outputed.
% 
% The table generated here is used in our paper:
%    "Active-set identification for interior point methods using controlled perturbations"
%
% Note: 
%     copy and paste the following line to line 73 and line 99 of pipm.m
%     to gother info we need:
%     
%      iter_info(p.counter.iterN+1,:) = [p.counter.iterN p.getMu p.prob.x']; save('iter_info.mat','iter_info');
%
clc;
clear all;

load verify.mat

param.mu_cap = 1e-06;
param.maxIter = 5;
param.verbose = 0;

% with per
% disp('With per')
param.iPer = 0.01;
p1 = pipm(A,b,c,param);
p1.solve;

load iter_info.mat;
iter_info_per = iter_info;

% without per
% disp('Without per')
param.iPer = 0;
p2 = pipm(A,b,c,param);
p2.solve;

load iter_info.mat;
iter_info_unp = iter_info;

title = {'Iter' 'MU' 'X1' 'X2'};

clearvars p1 p2 param iter_info; delete iter_info.mat;

% print table
fprintf('\\begin{table}\n\\label{tab-demo_lp}\n');
fprintf('\\caption{An example of predicting optimal active-set using perturbations}\n');
fprintf('\\begin{tabular}{c | c c >{\\columncolor{light-gray}}c c c | c c >{\\columncolor{light-gray}}c }\n\t\\hline\\hline\n');
fprintf('\t%4s & %9s & %9s & %9s & %14s & %14s & %9s & %9s & %9s \\\\ \n',...
    'Iter', '$\mu$', '$x_{1}$', '$x_{2}$', '$\lambda_{1}$', '$\lambda_{2}$', '$\mu$',...
    '$x_{1}$', '$x_{2}$');

fprintf('\t\\hline\n');

for i =1:6
    fprintf('\t%4d & %9.2e & %9.2e & %9.2e & %14.2e & %14.2e & %9.2e & %9.2e & %9.2e \\\\ \n',...
        iter_info_per(i,1), iter_info_per(i,2), iter_info_per(i,3), iter_info_per(i,4),iter_info_per(i,5),iter_info_per(i,6),...
        iter_info_unp(i,2), iter_info_unp(i,3), iter_info_unp(i,4));
end
fprintf('\t\\hline\\hline\n');
fprintf('\t\\end{tabular}\n');
fprintf('\\end{table}\n');