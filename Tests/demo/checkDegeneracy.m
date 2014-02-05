% This script is used to check if the test problems generate by the
% function generateRandomProb is primal nondegenerate.
clc;
clear;

numTestProb = 100;
i = 1;
primal_degen = zeros(numTestProb,1);

while i <= numTestProb
    
    [A, b, c] = generateRandomProb('m_min',10,'m_max',200,...
        'n_min',20,'n_max',500);
    
    [m, n]  = size(A);  A  = full(A);
    lb = zeros(n, 1);   ub = inf*ones(n, 1);
    
    % options = optimset('LargeScale', 'off', 'Algorithm','simplex','Display','off');
    options = optimset('Algorithm','active-set','Display','off');
    
    [x,Fval,Exitflag] = linprog(c,[],[],A,b,lb,ub,[],options);
    
    if Exitflag ~= 1
        continue;
    end
    
    if sum(x > 0) < m
        primal_degen(i) = 1;
    end
    
    clc;
    progress_bar_length = 100;
    current_bar_length = floor(progress_bar_length*i/numTestProb);
    output = [repmat('=',1,current_bar_length) '>' repmat(' ', 1, progress_bar_length-current_bar_length)];
    fprintf('[%s]  - %s %% done.\n',output,num2str(100*i/numTestProb));
    
    i = i+1;
    
end

disp('# of primal degenerate: ');
disp(num2str(sum(primal_degen)));