clc;
clear;

disp('[1] random test');
disp('[2] random degenerate')
randomTest = input('Choose test [1-2]: ');

disp(['Start the test ' num2str(randomTest)]);
%% Threshold
% C = [1e-06 1e-05 1e-04];
C = 1e-05;
K = [5 8 10 12 16];

% C = 1e-05
numTestProb = 5;
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
for i = 1:length(C)
    %% Reset the random number generator
    rng('default');
    rng(1);
    
    diff_ipm_splx = 0;
    diff_per_ipm = 0;
    diff_unper_splx = 0;
    diff_p_ip_sl = 0;
    
    rel_diff_ipm_splx = zeros(numTestProb,1);
    rel_diff_per_ipm = rel_diff_ipm_splx;
    rel_diff_per_splx = rel_diff_ipm_splx;
    
    fprintf('\n\nC = %9.2e\n\n', C(i));
    
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
        
        fprintf('%4s %4s %4s %2s %9s | %8s %10s %10s %10s | %10s %10s %10s %10s %13s %10s %13s\n',...
            'Prob', 'm', 'n', 'k', 'lambda',...
            'ac_p_ipm', 'ac_p_splx', 'ac_up_ipm', 'ac_up_splx',...
            'd_p_ip_sl',...
            'd_per_ipm', 'rd_per_ipm', ...
            'd_per_splx', 'rd_per_splx',...
            'd_ipm_splx', 'rd_ipm_splx');
        
        for k = 1:length(K)           
            %% per
            params_per.iPer = 1e-02;
            params_per.maxIter = K(k);
            per = pipm(A,b,c,params_per); per.solve;
            
            % build new perturbed problem
            current_p = per.perturbations.lambda * 1e-02;
            bp = b+A*current_p;
            cp = c+current_p;
            
            % solve this perturbed problem using ipm
            [m, n]  = size(A);   A  = full(A);
            lb = zeros(n, 1);  ub = inf*ones(n, 1);
            [p_ipm] = linprog(cp,[],[],A,bp,lb,ub,[],optn_ipm);
            
            ac_p_ipm = find(p_ipm < C(i));
            
            % solve this perturbed problem using simplex
            [p_splx] = linprog(cp,[],[],A,bp,lb,ub,[],optn_splx);
            
            ac_p_splx = find(p_splx < C(i));
            
            %% IPM solve for the PD
            [x_ipm] = linprog(c,[],[],A,b,lb,ub,[],optn_ipm);
            
            ac_up_ipm = find(x_ipm < C(i));
            
            %% Simplex
            [x_simplex,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog(c,[],[],A,b,lb,ub,[],optn_splx);
            ac_up_splx = find(x_simplex < C(i));
            
            %% stat
            diff_per_ipm = setdiff(ac_p_ipm, ac_up_ipm);
            diff_unper_splx = setdiff(ac_p_ipm, ac_up_splx);
            diff_ipm_splx = setdiff(ac_up_ipm, ac_up_splx);
            diff_p_ip_sl = setdiff(ac_p_ipm, ac_p_splx);
            
            chk_per_ipm  = isempty(diff_per_ipm);
            chk_per_splx = isempty(diff_unper_splx);
            chk_ipm_splx   = isempty(diff_ipm_splx);

            
            diff_per_ipm = union( diff_per_ipm, setdiff(ac_up_ipm, ac_p_ipm));
            diff_unper_splx = union( diff_unper_splx, setdiff(ac_up_splx, ac_p_ipm));
            diff_ipm_splx = union( diff_ipm_splx, setdiff(ac_up_splx, ac_up_ipm) );
            diff_p_ip_sl = union( diff_p_ip_sl, setdiff(ac_p_splx, ac_p_ipm ));
            
            rel_diff_per_ipm(j) = length(diff_per_ipm)/length(ac_up_ipm);
            rel_diff_per_splx(j) = length(diff_unper_splx)/length(ac_up_splx);
            rel_diff_ipm_splx(j) = length(diff_ipm_splx)/length(ac_up_splx);
            
            fprintf('%4d %4d %4d %2d %9.2e | %8d %10d %10d %10d | %10d %10d %10.3f %10d %13.3f %10d %13.3f\n',...
                j, m, n, K(k), mean(current_p),....
                length(ac_p_ipm), length(ac_p_splx), length(ac_up_ipm), length(ac_up_splx),...
                length(diff_p_ip_sl),...
                length(diff_per_ipm), rel_diff_per_ipm(j),...
                length(diff_unper_splx), rel_diff_per_splx(j),...
                length(diff_ipm_splx), rel_diff_ipm_splx(j));
            disp(chk_per_splx)
        end
        disp(' ')
        %         fprintf('--------------------------------------------------------------------------------------------------------------------------------------------\n')
        %         fprintf('Average %45s %13.3f %13s %13.3f %10s %13.3f\n', ' ',...
        %             mean(rel_diff_unper_ipm), ' ', mean(rel_diff_unper_splx),' ', mean(rel_diff_ipm_splx));
    end
end
diary off;