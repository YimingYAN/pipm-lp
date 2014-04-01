warning off;
clc;

list = readProbSet('netlib_IF_ID_CF.txt');
% list = readProbSet('netlib_crossover_IF_ID_CF.txt');



params.verbose          = 0;
params.iPer             = 0; % Unper pathfollow
params.doCrossOver      = 0;
params.mu_cap           = 1e-012;
params.tol              = 1e-015;

params_IF = params;
params_IF.actvPredStrtgy   = 'conservidfunc';

params_ID = params;
params_ID.actvPredStrtgy   = 'conservindica';

params_CF = params;
params_CF.actvPredStrtgy   = 'conservCutoff';


fprintf('          &    & M-8 & M-7 & M-6 & M-5 & M-4 & M-3 & M-2 & M-1 & M \\\\ \n');
disp('\hline')
for i = 1 : length(list)
%     load(list{i});
    
%     [A,b,c,FEASIBLE]=myPreprocess(A,b,c,lbounds,ubounds,BIG);
[A,b,c] = generateRandomProb;

NAME = ['        rdp_' int2str(i)];
    
%     [actualActv, exitflag, x] = solveLinprog(A,b,c, 'ipm');
    
    first_solve =  pipm(A,b,c, params); first_solve.solve;
    
    M = first_solve.getIPMIterCount; actualActv = find(first_solve.getx < 1e-05);
    
    str_IF = [NAME(end-8:end)   ' & IF']; 
    str_ID = [NAME(end-8:end)   ' & ID'];  
    str_CF = [' (M = ' int2str(M) ') & CF'];
    
    for k = M - 8 : 1 : M
        if k < 0
            str_ID = '--';
            str_CF = '--';
            continue;
        end
        
        params_IF.maxIter = k;
        params_ID.maxIter = k;
        params_CF.maxIter = k;
        
        p_IF = pipm(A,b,c, params_IF); p_IF.solve;
        p_ID = pipm(A,b,c, params_ID); p_ID.solve;
        p_CF = pipm(A,b,c, params_CF); p_CF.solve;
        
%         str_IF = [str_IF ' & ' sprintf('%3d', round(100*length(p_IF.getActv)/length(actualActv)))];
%         str_ID = [str_ID ' & ' sprintf('%3d', round(100*length(p_ID.getActv)/length(actualActv)))];
%         str_CF = [str_CF ' & ' sprintf('%3d', round(100*length(p_CF.getActv)/length(actualActv)))];
        
        [~, ~, cr_IF] = getCorrectionRatio(p_IF.getActv,actualActv);
        [~, ~, cr_ID] = getCorrectionRatio(p_ID.getActv,actualActv);
        [~, ~, cr_CF] = getCorrectionRatio(p_CF.getActv,actualActv);
        
        str_IF = [str_IF ' & ' sprintf('%3d', round(100*cr_IF))];
        str_ID = [str_ID ' & ' sprintf('%3d', round(100*cr_ID))];
        str_CF = [str_CF ' & ' sprintf('%3d', round(100*cr_CF))];
    end
    
    str_IF = [str_IF '\\ '];
    str_ID = [str_ID '\\ '];
    str_CF = [str_CF '\\ '];
    
    %disp(str_IF);
    disp(str_ID);
    disp(str_CF);
    disp('\hline')
    
    
    
end