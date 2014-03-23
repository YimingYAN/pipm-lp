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

[x,fval,exitflag] = fmincon('fun_gamma', t0, [], [], Aeq, beq,lb, ub, [], options);

gamma = fval
