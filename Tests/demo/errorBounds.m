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