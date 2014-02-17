% This script is used to illustarte that we could find a vector of
% perturbations such that epsilon(A,b,c) < epsilon(A,b_hat,c_hat)
% 
% Used for demostration.
%
% Yiming

clc;
clear;

load verify;

% lambda = [0.01; 0.05];
lambda = rand*[1; 5];
disp('    lambda:') 
disp(lambda)
disp(' ')

disp('    optx   opts')
disp([optx opts]);

a = optx < 1e-05;
i = optx > 1e-05;

A_a = A( :, a ); 
A_i = A(:, i); 

chk1 = lambda( i ) + inv(A_i)*A_a*lambda( a ) > 0;
chk2 = lambda( a ) - (inv(A_i)*A_a)'*lambda( i ) > 0;

disp('lambda( i ) + inv(A_i)*A_a*lambda( a ) > 0 ? ')
if chk1 disp('Yes'); else disp('No'); end; disp(' ')
disp('lambda( a ) - transpose((inv(A_i)*A_a))*lambda( i ) > 0 ?')
if chk2 disp('Yes'); else disp('No'); end; disp(' ')

epsilon = min(min(optx(i)),min(opts(a)));


b_l = b + A*lambda;
c_l = c+lambda;

lb = zeros(length(optx),1);

[p,~,~,~,dual] = linprog(c_l,[],[],A,b_l,lb);
q = dual.lower;
disp('    p         q')
disp([p q]);

epsilon_hat = min(min(p(i)),min(q(a)));

disp('    epsilon   epsilon_hat ')
disp([epsilon epsilon_hat]);