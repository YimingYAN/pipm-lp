function z = dist_x(t)
% t = (x1,x2,s1,s2,y)
% z is a real number
% x* = (1,0)

x_star = [1; 0];
w = [-t(1:4); t(1:2)'*t(3:4)];

z = -norm(t(1:2) - x_star) / ( norm( min(t(1:2), t(3:4)) ) + norm( w(w>0) ) );