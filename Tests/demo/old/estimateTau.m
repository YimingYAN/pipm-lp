% This script is used to estimate the constant tau_p needed for
% the error bounds.
% Random problem is tested.

% Uncomment iter_x(:,p.counter.iterN+1) = p.prob.x; iter_s(:,p.counter.iterN+1) = p.prob.s; iter_y(:,p.counter.iterN+1) = p.prob.y; save('iter_xys.mat', 'iter_x', 'iter_y', 'iter_s');
% in pipm

function estimateTau()

clc;
close all;

load verify6.mat

param.mu_cap = 1e-09;
param.tol = 1e-09;
param.verbose = 0;
param.doCrossOver = 0;

% with per
% disp('With per')
param.iPer = 1e-02;
p1 = pipm(A,b,c,param);
p1.solve;

load iter_xys.mat;
X = iter_x; Y = iter_y; S = iter_s;

% X: each column of X contians iterate of x

%% estimate tau_p
L = size(X,2);
tau_p = zeros(L,1);

for i =1: L
    eb = errorBounds(A, b, c, X(:,i), Y(:,i), S(:,i));
    
    tau_p(i) = norm( X(:,i) - optx)/eb;
end

tau_p

tau_p = mean(tau_p)

delete iter_xys.mat

end