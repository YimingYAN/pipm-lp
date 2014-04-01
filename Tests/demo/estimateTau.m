% find the value of tau_p and tau_d numerically for TEST PROBLEM
% verify.
%
% Essentially we are solving
%           max ||x-x*||/(r(x,s) + w(x,s))
%           s.t. Ax = b, A'y+s = c, x + lambda > 0, s + lambda > 0
% and 
%           max ||s-s*||/(r(x,s) + w(x,s))
%           s.t. Ax = b, A'y+s = c, x + lambda > 0, s + lambda > 0
 
% For fmincon:
%           dist_x: -norm(x-x*)/(r(x,s) + w(x,s))
%           dist_y: -norm(s-s*)/(r(x,s) + w(x,s))
%                t: (x1,x2,s1,s2,y)
%
clear;
clc;
close all;


load verify;

[m,n] = size(A);

Aeq = [ A zeros(m,n) zeros(m,m);
    zeros(n,n) eye(n) A'];

beq = [b;c];

t0 = [1;1;1;1;1];

ub = inf*ones(5,1);

options = optimset('Display', 'off','Algorithm','active-set');

lambda = 1e-02;

lb = [-lambda; -5*lambda; -lambda; -5*lambda; -inf];

[x,fval1,exitflag] = fmincon('dist_x', t0, [], [], Aeq, beq,lb, ub, [], options);
[x,fval2,exitflag] = fmincon('dist_s', t0, [], [], Aeq, beq,lb, ub, [], options);

tau_p = -fval1;
tau_d = -fval2;

disp('tau_p :'); disp(tau_p)
disp('tau_d :'); disp(tau_d)