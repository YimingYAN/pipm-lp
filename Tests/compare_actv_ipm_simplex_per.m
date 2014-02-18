clc;
clear;

disp('[1] random test');
disp('[2] random degenerate')
randomTest = input('Choose test [1-2]: ');
%% Threshold
C = [1e-06 1e-05 1e-04 1e-03];

C = 1e-05;
numTestProb =20;
!rm compare_actv_ipm_simplex_per_demo.log
diary('compare_actv_ipm_simplex_per_demo.log');

%% With perturbations
params_per.verbose        = 0;
params_per.iPer           = 1e-02;
params_per.actvPredStrtgy = 'conservCutoff';
params_per.doCrossOver    = 0;
params_per.mu_cap         = 1e-032;    % Avoid being terminated by mu_cap
params_per.tol            = 1e-032;

for i = 1:length(C)
    fprintf('\n\n========== C = %9.2e ==========\n\n', C(i));
    
    for k=12:18
        params_per.maxIter = k;
        fprintf('\n\nk = %d\n\n', k);
        %% Reset the random number generator
        rng('default');
        rng(1);
        
        diff_ipm_splx = 0;
        diff_per_ipm = 0;
        diff_per_splx = 0;
        
        rel_diff_ipm_splx = zeros(numTestProb,1);
        rel_diff_per_ipm = rel_diff_ipm_splx;
        rel_diff_per_splx = rel_diff_ipm_splx;
        
        mu = rel_diff_ipm_splx;
        
        
        
        fprintf('%4s %4s %9s %10s %10s %10s %10s %13s %13s %13s %10s %13s %10s %10s %10s\n',...
            'Prob', 'n', 'mu', 'actv_per', 'actv_ipm', 'actv_splx',...
            'd_per_ipm', 'r_d_per_ipm', ...
            'd_per_splx', 'r_d_per_splx',...
            'd_ipm_splx', 'r_d_ipm_splx',...
            'per<ipm', 'per<splx', 'ipm<splx');
        
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
            %% Unper
            per = pipm(A,b,c,params_per); per.solve;
            x_per = per.prob.x;
            actv_per = find(x_per < C(i));
            
            
            %% IPM
            n  = size(A, 2);   A  = full(A);
            lb = zeros(n, 1);  ub = inf*ones(n, 1);
            
            options=[];
            options = optimset('Algorithm', 'interior-point','Display','off');
            
            [x_ipm] = linprog(c,[],[],A,b,lb,ub,[],options);
            
            actv_ipm = find(x_ipm < C(i));
            
            %% Simplex
            options = [];
            options = optimset('LargeScale', 'off', 'Algorithm','simplex','Display','off');
            %         options = optimset('Algorithm','active-set','Display','off');
            
            [x_simplex,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog(c,[],[],A,b,lb,ub,[],options);
            
            %        actv_simplex = find(x_simplex == 0);
            %         actv_simplex = find(x_simplex < 1e-05);
            actv_splx = find(x_simplex < C(i));
            
            %% stat
            diff_per_ipm = setdiff(actv_per, actv_ipm);
            diff_per_splx = setdiff(actv_per, actv_splx);
            diff_ipm_splx = setdiff(actv_ipm, actv_splx);
            
            chk_per_ipm  = isempty(diff_per_ipm);
            chk_per_splx = isempty(diff_per_splx);
            chk_ipm_splx   = isempty(diff_ipm_splx);
            
            diff_per_ipm = union( diff_per_ipm, setdiff(actv_ipm, actv_per));
            diff_per_splx = union( diff_per_splx, setdiff(actv_splx, actv_per));
            diff_ipm_splx = union( diff_ipm_splx, setdiff(actv_splx, actv_ipm) );
            
            rel_diff_per_ipm(j) = length(diff_per_ipm)/length(actv_ipm);
            rel_diff_per_splx(j) = length(diff_per_splx)/length(actv_splx);
            rel_diff_ipm_splx(j) = length(diff_ipm_splx)/length(actv_splx);
            
            mu(j) = per.getMu;
            
            fprintf('%4d %4d %9.2e %10d %10d %10d %10d %13.3f %13d %13.3f %10d %13.3f %10d %10d %10d\n',...
                j, n, per.getMu, length(actv_per), length(actv_ipm), length(actv_splx),...
                length(diff_per_ipm), rel_diff_per_ipm(j),...
                length(diff_per_splx), rel_diff_per_splx(j),...
                length(diff_ipm_splx), rel_diff_ipm_splx(j),...
                chk_per_ipm, chk_per_splx, chk_ipm_splx);
        end
        fprintf('--------------------------------------------------------------------------------------------------------------------------------------------\n')
        fprintf('Average   %9.2e %43s %13.3f %13s %13.3f %10s %13.3f\n', ...
            mean(mu),' ', ...
            mean(rel_diff_per_ipm), ' ', mean(rel_diff_per_splx),' ', mean(rel_diff_ipm_splx));
    end
end
diary off;