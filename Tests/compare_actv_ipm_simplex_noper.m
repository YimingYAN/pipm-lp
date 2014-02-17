clc;
clear;

disp('[1] random test');
disp('[2] random degenerate')
randomTest = input('Choose test [1-2]: ');
%% Threshold
C = [1e-06 1e-05 1e-04 1e-03];
% C = 1e-05
numTestProb =20;
!rm compare_actv_ipm_simplex_noper_demo.log
diary('compare_actv_ipm_simplex_noper_demo.log');

%% Without perturbations
params_unper.verbose        = 0;
params_unper.iPer           = 0;
params_unper.actvPredStrtgy = 'conservCutoff';
params_unper.doCrossOver    = 0;
params_unper.mu_cap         = 1e-09;
params_unper.tol            = 1e-09;

for i = 1:length(C)
    %% Reset the random number generator
    rng('default');
    rng(1);
    
    diff_ipm_splx = 0;
    diff_unper_ipm = 0;
    diff_unper_splx = 0;
    
    rel_diff_ipm_splx = zeros(numTestProb,1);
    rel_diff_unper_ipm = rel_diff_ipm_splx;
    rel_diff_unper_splx = rel_diff_ipm_splx;
    
    fprintf('\n\nC = %9.2e\n\n', C(i));
    fprintf('%4s %4s %10s %10s %10s %10s %13s %13s %13s %10s %13s %10s %10s %10s\n',...
        'Prob', 'n', 'actv_unper', 'actv_ipm', 'actv_splx',...
        'd_unp_ipm', 'r_d_unp_ipm', ...
        'd_unp_splx', 'r_d_unp_splx',...
        'd_ipm_splx', 'r_d_ipm_splx',...
        'unper<ipm', 'unper<splx', 'ipm<splx');
    
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
        unper = pipm(A,b,c,params_unper); unper.solve;
        x_unper = unper.prob.x;
        actv_unper = find(x_unper < C(i));
        
        
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
        diff_unper_ipm = setdiff(actv_unper, actv_ipm);
        diff_unper_splx = setdiff(actv_unper, actv_splx);
        diff_ipm_splx = setdiff(actv_ipm, actv_splx);
        
        chk_unper_ipm  = isempty(diff_unper_ipm);
        chk_unper_splx = isempty(diff_unper_splx);
        chk_ipm_splx   = isempty(diff_ipm_splx);
        
        diff_unper_ipm = union( diff_unper_ipm, setdiff(actv_ipm, actv_unper));
        diff_unper_splx = union( diff_unper_splx, setdiff(actv_splx, actv_unper));
        diff_ipm_splx = union( diff_ipm_splx, setdiff(actv_splx, actv_ipm) );
        
        rel_diff_unper_ipm(j) = length(diff_unper_ipm)/length(actv_ipm);
        rel_diff_unper_splx(j) = length(diff_unper_splx)/length(actv_splx);
        rel_diff_ipm_splx(j) = length(diff_ipm_splx)/length(actv_splx);
        
        fprintf('%4d %4d %10d %10d %10d %10d %13.3f %13d %13.3f %10d %13.3f %10d %10d %10d\n',...
            j, n, length(actv_unper), length(actv_ipm), length(actv_splx),...
            length(diff_unper_ipm), rel_diff_unper_ipm(j),...
            length(diff_unper_splx), rel_diff_unper_splx(j),...
            length(diff_ipm_splx), rel_diff_ipm_splx(j),...
            chk_unper_ipm, chk_unper_splx, chk_ipm_splx);
    end
    fprintf('--------------------------------------------------------------------------------------------------------------------------------------------\n')
    fprintf('Average %45s %13.3f %13s %13.3f %10s %13.3f\n', ' ',...
        mean(rel_diff_unper_ipm), ' ', mean(rel_diff_unper_splx),' ', mean(rel_diff_ipm_splx));
end
diary off;