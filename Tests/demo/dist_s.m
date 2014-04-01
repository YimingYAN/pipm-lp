function z = dist_s(t)
% t = (x1,x2,s1,s2,y)
% z is a real number
% s* = (0,1)

s_star = [0; 1];
w = [-t(1:4); t(1:2)'*t(3:4)];

z = -norm(t(3:4) - s_star) / ( norm( min(t(1:2), t(3:4)) ) + norm( w(w>0) ) );