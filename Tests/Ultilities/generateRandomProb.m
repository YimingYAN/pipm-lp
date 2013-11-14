function [A, b, c] = generateRandomProb(varargin)
% m --- (10,200)
% n --- (20,500)
% density --- (0.4,0.8)

debug = 0;
m_min = 10; m_max = 200;
n_min = 20; n_max = 500;

for i=1:nargin
    switch lower(varargin{i})
        case {'debug'}
            debug = 1;
        case {'m_min'}
            m_min = varargin{i+1};
            i=i+1;
        case {'m_max'}
            m_max = varargin{i+1};
            i = i+1;
        case {'n_min'}
            n_min = varargin{i+1};
            i = i+1;
        case {'n_max'}
            n_max = varargin{i+1};
            i = i+1;
    end
end

if n_max - n_min <=0 || m_max-m_min <=0
    error('make sure n_max > n_min and m_max > m_min');
end

m = zeros(1, 1, 'int32');
n = zeros(1, 1, 'int32');
while m >= n || n <= 2*m || n >= 7*m
    % get m
    m = rand;
    m = m_min+floor((m_max-m_min)*m);
    
    % get n
    n = rand;
    n = n_min+floor((n_max-n_min)*n);
end
% get density
density = rand;
density = 0.4+0.4*density;

% generate a sparse random matrix of given density
A = sprandn(m,n,density);

% choose feasible x, y, s at random, with x and s each about half-full
xfeas = [rand(floor(n/2),1); zeros(n-floor(n/2),1)];
sfeas = [zeros(floor(n/2),1); rand(n-floor(n/2),1)];
xfeas = xfeas(randperm(n)); sfeas = sfeas(randperm(n)); 
yfeas = (rand(m,1)-0.5)*4;

% choose b and c to make this (x,y,s) feasible
b = A*xfeas; c=A'*yfeas+sfeas;
end