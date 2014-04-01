% This scripts tries to find out the largest value of lambda that will
% maintain the optimal actv
% for test problem verify.mat

clear;
clc;

load verify;

[m,n] = size(A);

lb = zeros(n,1);



range = -10:10; range = 10.^range;

for i = 1 : length(range)
    lambda = range(i);
    
    bp = b + A*[lambda; 5*lambda];
    cp = c + [lambda; 5*lambda];
    
    [actualActv, exitflag, x] = solveLinprog(A,bp,cp, 'splx');
x    
end
