clear all;
clc;
close all;

numTestProb =3;


param.iPer = 0;
param.tol = 1e-9;
param.mu_cap = 1e-9;
param.verbose = 0;
param.doCrossOver = 0;

rng(1);

for i = 1 : numTestProb
    % random problems
    [A,b,c] = generateRandomProb;

    
    % with per
    % disp('With per')
    up = pipm(A,b,c,param);
    up.solve;
    
    optx = up.getx;
    
    param.iPer = 1e-02;
    p = pipm(A,b,c,param);
    p.solve;
    
    load iter_xys.mat;
    X = iter_x; Y = iter_y; S = iter_s;
    
    % X: each column of X contians iterate of x
    
    %% estimate tau_p
    L = size(X,2);
    tau_p = zeros(L,1);
    eb = zeros(L,1);
    
    for k =1: L
        eb(k) = errorBounds(A, b, c, X(:,k), Y(:,k), S(:,k));
        tau_p(k) = norm( X(:,k) - optx)/eb(k);
        
    end
    
    
    
    [tau_p eb]
    tau_p = mean(tau_p)
    
    delete iter_xys.mat
    
end
