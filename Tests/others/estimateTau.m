% This script is used to estimate the constant tau_p needed for
% the error bounds.
% Random problem is tested.

% Uncomment iter_x(:,p.counter.iterN+1) = p.prob.x; iter_s(:,p.counter.iterN+1) = p.prob.s; iter_y(:,p.counter.iterN+1) = p.prob.y; save('iter_xys.mat', 'iter_x', 'iter_y', 'iter_s');
% in pipm

function estimateTau()

clc;
close all;

load verify.mat

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

function eb = errorBounds(A, b, c, x, y, s)

% We assume the LP has a unique solution.
%
% We do not ristrict the points (xk,yk,sk) to satisfy
%     Ax=b, A'y+s=c, x>=0, s>=0
%
% Version 0.1
% Author: Yiming Yan


% initialize ------------
[m,n] = size(A);
res_p = A*x-b;
res_d = c-A'*y;

% get y- and y+ ------------
ymns = -y; yplus = zeros(m,1);
indx = y >= 0;
yplus(indx) = y(indx); ymns(indx) = 0;

% get error bounds ------------
% r
r_1 = min(x, res_d);
r_2 = min(yplus, res_p);
r_3 = min(ymns, -res_p);
r = [r_1; r_2; r_3];
r = norm(r,2);

% w
w = [-res_d; -res_p; res_p; -x; x'*res_d+y'*res_p];
indx2 = w < 0;
w(indx2) = [];
w = norm(w,2);

eb = r+w;
end