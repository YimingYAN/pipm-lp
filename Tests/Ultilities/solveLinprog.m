function [actualActv, exitflag, x] = solveLinprog(A,b,c, alg)
%SOLVELINPROG Make use of the MATLAB's linear programming solver: LINPOG
%   alg: 
%       splx     -  simplex solver
%       ipm      -  interior-point solver
n  = size(A, 2);   A  = full(A);
lb = zeros(n, 1);  ub = inf*ones(n, 1);

switch lower(alg)
    case 'splx'
        options = optimoptions('linprog',...
            'Algorithm','dual-simplex',...
            'TolFun', 1e-09,...
            'Display','off');
    case 'ipm'
        options = optimoptions('linprog',...
            'Algorithm','interior-point',...
            'TolFun', 1e-09,...
            'Display','off');
    otherwise
        error('solveLinprog: no such a solver');
end

[xsol,~,exitflag, output] = linprog(c,[],[],A,b,lb,ub,[],options);

% Get actual actv
actualActv = find(xsol<1e-05);

if nargin > 2
    x = xsol;
end
end
