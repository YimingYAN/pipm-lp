function compare_original_sol_interior_point_actv_set()
% compare_original_sol_interior_point_actv_set This function compares
% difference between the optimal solution obtained from interior point
% solver and that from the active-set solution for random and ramdom
% degenerate problems. We solve the original problem here.
%
% March 28 2014
% Yiming Yan

numTestProb = 100;

rng('default');
rng(1);

str{1} = 'Random       :[';
str{2} = 'Random Degen :[';

title_length = length(str{1});
current_position = title_length + 2;

str{1} = ['Random       : [' repmat(' ', 1, numTestProb) ']'];
str{2} = ['Random Degen : [' repmat(' ', 1, numTestProb) ']'];
actv_diff_random = zeros(numTestProb,1);
actv_diff_rnddegen = zeros(numTestProb,1);

for i = 1:numTestProb
    str{1}(current_position) = '-';
    clc; disp(str{1}); 
    
    [A,b,c] = generateRandomProb;
    
    [actualActv_ipm, exitflag, x] = solveLinprog(A,b,c, 'ipm');
    [actualActv_actv, exitflag, x] = solveLinprog(A,b,c, 'actv-set');
    
    str{1}(current_position) = '=';
    clc; disp(str{1}); 
    
    tmp =  union( setdiff(actualActv_ipm, actualActv_actv),  setdiff(actualActv_actv, actualActv_ipm) );
    actv_diff_random(i) = length(tmp)/length(actualActv_actv);

    current_position = current_position +1;
end

current_position = title_length + 2;
rng('default');
rng(1);

for i=1:numTestProb
    clc; disp(str{1});
    str{2}(current_position) = '-';
    disp(str{2});
    
    [A,b,c] = generateDegenProb;
    [actualActv_ipm, exitflag, x] = solveLinprog(A,b,c, 'ipm');
    [actualActv_actv, exitflag, x] = solveLinprog(A,b,c,'actv-set');
    
    tmp =  union( setdiff(actualActv_ipm, actualActv_actv),  setdiff(actualActv_actv, actualActv_ipm) );
    actv_diff_rnddegen(i) = length(tmp)/length(actualActv_actv);
    
    clc; disp(str{1}); 
    str{2}(current_position) = '=';
    disp(str{2})
    
    current_position = current_position +1;
end

fprintf('%10s %10s %15s\n', ' ', 'random', 'random_degen' );
fprintf('%10s %10.2f %10.2f\n', 'Avg_diff', mean(actv_diff_random), mean(actv_diff_rnddegen) );

end