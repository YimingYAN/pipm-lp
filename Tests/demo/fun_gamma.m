function z = fun_gamma(t)

lambda = [1e-02; 5*1e-02];

z = 2*min((t(1:2) + lambda).*(t(3:4) + lambda)) / (t(1:2)'*t(3:4));


end