% find the value of tau_p and tau_d numerically for problem
% verify.
%
% Essentially we are solving
%           max ||x-x*||/(r(x,s) + w(x,s))
%           s.t. Ax = b, A'y+s = c, x + lambda > 0, s + lambda > 0
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
range = -10:10; range = 10.^range;

for i = 1 : length(range)
    lambda = range(i);
 
    lb = [-lambda; -5*lambda; -lambda; -5*lambda; -inf];
    
    [x,fval1,exitflag] = fmincon('dist_x', t0, [], [], Aeq, beq,lb, ub, [], options);
    [x,fval2,exitflag] = fmincon('dist_s', t0, [], [], Aeq, beq,lb, ub, [], options);
    
    tau_p(i) = -fval1;
    tau_d(i) = -fval2;
end

figure;
fig_x = plot(log10(range),tau_p,'-*b');
xlabel('log10($\lambda$)','interpreter','latex');
ylabel('$\tau_{p}$','interpreter','latex');
title('$\tau_{p}$ versus $\lambda$','interpreter','latex')

figure;
fig_s = plot(log10(range),tau_d,'-*r');
xlabel('log10($\lambda$)','interpreter','latex');
ylabel('$\tau_{d}$','interpreter','latex');
title('$\tau_{d}$ versus $\lambda$','interpreter','latex')