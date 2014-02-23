% This script is used to verify the difference bwtween the actual
% acrive-set of the perturbed porblem and the actual active-set of the
% original problem. We are trying to answer the following questions:
%
% 1. Are the actual active sets the same?
%
% 2. Do the pertuebed problems have unique and nodegenerate solution?
%
% 3. Is it true when the perturbation are small enough perturbed problems
% have the same active-set as the original problem?
%
% Date  : 23 Feb 2014
% Author: Yiming Yan

clc;
clear;

disp('[1] random test');
disp('[2] random degenerate')
randomTest = input('Choose test [1-2]: ');

disp(['Start the test ' num2str(randomTest)]);
%% Threshold
% C = [1e-06 1e-05 1e-04];
C = 1e-05;
K = [0 10 18];    
% the second number is chosen by having largest gap between pac and oac 
                  

% C = 1e-05
numTestProb = 100;

!rm check_ac_p_ipm_vs_actv_unper.log
diary('check_ac_p_ipm_vs_actv_unper.log');

%% Without perturbations
params_per.verbose        = 0;
params_per.actvPredStrtgy = 'conservCutoff';
params_per.doCrossOver    = 0;
params_per.mu_cap         = 1e-09;
params_per.tol            = 1e-09;

%% linprog
optn_ipm = optimset('Algorithm', 'interior-point','Display','off');
optn_splx = optimset('LargeScale', 'off', 'Algorithm','simplex','Display','off');

%% main loop
%for i = 1:length(C)

%fprintf('\n\nC = %9.2e\n\n', C(i));

num_pac_ne_oac               = zeros(length(K) + 1, 1);
num_per_unique_nondeg        = num_pac_ne_oac;

rel_diff_up_ip_sl   = zeros(numTestProb, length(K) + 1);
rel_diff_p_up_ipm   = rel_diff_up_ip_sl;
rel_diff_p_up_splx  = rel_diff_up_ip_sl;
rel_diff_p_ip_sl    = rel_diff_up_ip_sl;
    

for k = 1:length(K) + 1
    %% Reset the random number generator
    rng('default');
    rng(1);
    
    diff_up_ip_sl   = 0;
    diff_p_up_ipm   = 0;
    diff_p_up_splx  = 0;
    diff_p_ip_sl    = 0;
    
    disp(' ')
    %fprintf('%4s %4s %4s %2s %9s | %8s %10s %10s %10s | %10s( %10s | %10s %10s | %10s %10s | %10s %10s\n',...
    fprintf('%4s %4s %4s %2s %9s | %8s %10s %10s %10s | %10s (%5s%%) %10s (%5s%%) %10s (%5s%%) %10s (%5s%%)\n',...
        'Prob', 'm', 'n', 'k', 'lambda',...
        'ac_p_ipm',     'ac_p_splx', 'ac_up_ipm', 'ac_up_splx',...
        'd_p_ip_sl',    ' ',...
        'd_up_ip_sl',   ' ',...
        'd_p_up_ip',    ' ', ...
        'd_p_up_sp',    ' ');
    
    for j=1:numTestProb
        if randomTest == 1
            [A, b, c] = ...
                generateRandomProb('m_min',10,'m_max',200,...
                'n_min',20,'n_max',500);
            
        elseif randomTest == 2
            
            [A, b, c] =...
                generateDegenProb('m_min',10,'m_max',200,...
                'n_min',20,'n_max',500);
        end
        
        [m, n]  = size(A);
        
        
        %% Solve the perturbed problem using linprog
        % get iper
        if k <= length(K)
            kth = K(k);
            params_per.iPer = 1e-02;
            params_per.maxIter = K(k);
            per = pipm(A,b,c,params_per); per.solve;
            current_p = per.perturbations.lambda;
            
        else
            current_p = 0*ones(n,1);
            kth = -1;
        end
        
        % build new perturbed problem
        
        bp = b + A*current_p;
        cp = c + current_p;
        A  = full(A);
        lb = zeros(n, 1);  ub = inf*ones(n, 1);
        
        % solve this perturbed problem using ipm
        
        [p_ipm] = linprog(cp, [], [], A, bp, lb, ub, [], optn_ipm);
        ac_p_ipm = find(p_ipm < C);
        
        % solve this perturbed problem using simplex
        [p_splx] = linprog(cp, [], [], A, bp, lb, ub, [], optn_splx);
        ac_p_splx = find( p_splx < C  );
        
        %% Solve the original problem using linprog
        % solve the original problme using IPM
        [x_ipm] = linprog(c,[],[],A,b,lb,ub,[],optn_ipm);
        ac_up_ipm = find(x_ipm < C );
        
        %% solve the original problme using simplex
        [x_simplex,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog(c,[],[],A,b,lb,ub,[],optn_splx);
        ac_up_splx = find(x_simplex < C );
        
        %% stat
        diff_p_ip_sl    = setdiff(ac_p_ipm, ac_p_splx);
        diff_up_ip_sl   = setdiff(ac_up_ipm, ac_up_splx);
        diff_p_up_ipm   = setdiff(ac_p_ipm, ac_up_ipm);
        diff_p_up_splx  = setdiff(ac_p_splx, ac_up_splx);
        
        chk_p_ip_sl     = isempty(diff_p_ip_sl);
        chk_up_ip_sl    = isempty(diff_up_ip_sl);
        chk_p_up_ipm    = isempty(diff_p_up_ipm);
        chk_p_up_splx   = isempty(diff_p_up_splx);
        
        diff_p_ip_sl    = union( diff_p_ip_sl,   setdiff(ac_p_splx,  ac_p_ipm ));
        diff_up_ip_sl   = union( diff_up_ip_sl,  setdiff(ac_up_splx, ac_up_ipm));
        diff_p_up_ipm   = union( diff_p_up_ipm,  setdiff(ac_up_ipm,  ac_p_ipm) );
        diff_p_up_splx  = union( diff_p_up_splx, setdiff(ac_up_splx, ac_p_splx));
        
        rel_diff_p_ip_sl(j,k)     = length(diff_p_ip_sl)  / length(ac_p_splx);
        rel_diff_up_ip_sl(j,k)    = length(diff_up_ip_sl) / length(ac_up_splx);
        rel_diff_p_up_ipm(j,k)    = length(diff_p_up_ipm) / length(ac_up_ipm);
        rel_diff_p_up_splx(j,k)   = length(diff_p_up_splx)/ length(ac_up_splx);
        
        if isempty(diff_p_ip_sl)
            num_per_unique_nondeg(k) = num_per_unique_nondeg(k) + 1;
        end
        
        if ~isempty(diff_p_up_ipm) || ~isempty(diff_p_up_splx)
            num_pac_ne_oac(k) = num_pac_ne_oac(k) + 1;
        end
        
        
        
        %% Output
        
        %fprintf('%4d %4d %4d %2d %9.2e | %8d %10d %10d %10d | %10d %10.3f | %10d %10.3f | %10d %10.3f | %10d %10.3f\n',...
        fprintf('%4d %4d %4d %2d %9.2e | %8d %10d %10d %10d | %10d (%5.2f%%) %10d (%5.2f%%) %10d (%5.2f%%) %10d (%5.2f%%)\n',...
            j, m, n, kth, mean(current_p),....
            length(ac_p_ipm),       length(ac_p_splx), length(ac_up_ipm), length(ac_up_splx),...
            length(diff_p_ip_sl),   100*rel_diff_p_ip_sl(j,k),...
            length(diff_up_ip_sl),  100*rel_diff_up_ip_sl(j,k),...
            length(diff_p_up_ipm),  100*rel_diff_p_up_ipm(j,k),...
            length(diff_p_up_splx), 100*rel_diff_p_up_splx(j,k) );
        %             disp(chk_per_splx)
    end
    fprintf('--------------------------------------------------------------------------------------------------------------------------------------------\n');
    fprintf('num_per_unique_nondeg (%%) = %d (%5.2f%%) \n', num_per_unique_nondeg(k), num_per_unique_nondeg(k)*100/numTestProb);
    fprintf('num_pac_ne_oac (%%)        = %d (%5.2f%%) \n', num_pac_ne_oac(k), num_pac_ne_oac(k)*100/numTestProb );
    fprintf('num_pac_eql_oac (%%)       = %d (%5.2f%%) \n', numTestProb - num_pac_ne_oac(k), 100*(1- num_pac_ne_oac(k)/numTestProb) );
end
%end

%% Summary table
disp('============')
fprintf('\n%9s %15s %15s %15s %15s %15s\n',....
    'lambda', '%_p_uni_ndeg', '%_pac_ne_oac',...
    'avg_r_p_ip_sl', 'avg_r_p_up_ip', 'avg_r_p_up_sp');
for i = 1: length(num_pac_ne_oac)
    fprintf('%9.2e %14.2f%% %14.2f%% %14.2f%% %14.2f%% %14.2f%%\n',...
        0, 100*num_per_unique_nondeg(i)/numTestProb, 100*num_pac_ne_oac(i)/numTestProb,...
        100*mean(rel_diff_p_ip_sl(:,i)),...
        100*mean(rel_diff_p_up_ipm(:,i)),...
        100*mean(rel_diff_p_up_splx(:,i)));
    
end

diary off;