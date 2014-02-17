clc;
clear;

%% Threshold
demo = 1;

% demo = 0;


%% setup
C = [1e-06 1e-05 1e-04 1e-03];
if demo
    numTestProb = 20;
    diary('compare_actv_ipm_simplex_noper_demo.log');
else
    numTestProb = 100;
    diary('compare_actv_ipm_simplex_noper.log');
    fprintf('\n%9s %9s\n', 'C', 'rel_diff_avg');
end
for i = 1:length(C)
    %% Reset the random number generator
    rng('default');
    rng(1);
    diff = 0;
    rel_diff = zeros(numTestProb,1);
    
    if demo
        fprintf('\n\nC = %9.2e\n\n', C(i));
        fprintf('%4s %4s %10s %12s %10s %10s %s\n',...
            'Prob', 'n', 'actv_ipm', 'actv_simplex', 'dff', 'rel_diff', 'ipm in simplex?');
    end
    
    for j=1:numTestProb
        [A, b, c] = ...
            generateRandomProb('m_min',10,'m_max',200,...
            'n_min',20,'n_max',500);
        
%         [A, b, c] =...
%             generateDegenProb('m_min',10,'m_max',200,...
%             'n_min',20,'n_max',500);
        
        n  = size(A, 2);   A  = full(A);
        lb = zeros(n, 1);  ub = inf*ones(n, 1);
        
        %% IPM
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
        actv_simplex = find(x_simplex < C(i));
        
        %% stat
        diff = setdiff(actv_ipm, actv_simplex);
        
        chk = isempty(diff);
        
        diff = union( diff, setdiff(actv_simplex, actv_ipm) );
        
        rel_diff(j) = length(diff)/length(actv_simplex);
        
        if demo
            fprintf('%4d %4d %10d %12d %10d %10.2f %10d\n',...
                j, n, length(actv_ipm), length(actv_simplex),...
                length(diff), rel_diff(j),chk );
        end
    end
    if demo
        fprintf('--------------------------------------------------------\n')
        fprintf('Average %36s %10.2f\n', ' ', mean(rel_diff));
    else
        fprintf('%9.2e %9.2f\n',C(i), mean(rel_diff));
    end
end
diary off;