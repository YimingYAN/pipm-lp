function [actualActv, exitflag, x] = solveLinprog(A,b,c, alg)
%SOLVELINPROG Make use of the MATLAB's linear programming solver: LINPOG
%   alg: 
%       splx     -  simplex solver
%       ipm      -  interior-point solver
%       actv-set -  active-set solver
n  = size(A, 2);   A  = full(A);
lb = zeros(n, 1);  ub = inf*ones(n, 1);

switch lower(alg)
    case 'splx'
        options = optimset('LargeScale', 'off', 'Algorithm','simplex','Display','off');
    case 'ipm'
        options = optimset('Algorithm', 'interior-point','Display','off');
    case 'actv-set'
        options = optimset('Algorithm','active-set','Display','off');
    otherwise
        error('solveLinprog: no such a solver');
end

[xsol,~,exitflag] = linprog(c,[],[],A,b,lb,ub,[],options);

% Get actual actv
actualActv = find(xsol<1e-05);

if nargin > 2
    x = xsol;
end
end
