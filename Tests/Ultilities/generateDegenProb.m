function [A, b, c] = generateDegenProb(varargin)
% m --- (m_min,m_max)
% n --- (n_min,n_max)
% density --- (0.4,0.8)
% degenDegree --- (0.1 0.3) 


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

% get degenerate degree
degenDegree = rand;
degenDegree = 0.1+0.2*degenDegree;

% generate a sparse random matrix of given density
A = sprandn(m,n,density); 

% get a basis
[q,r,e] = qr(A');
[junk,ind]=max(e);
ind = ind(1:size(A',2));
B = A(:,ind); 
N = A(:,setdiff(1:n,ind));
A = [B N];

% generate a degenerate optimal solution
x = [rand(floor(m*(1-degenDegree)),1); zeros(n-floor(m*(1-degenDegree)),1)];
y = (rand(m,1)-0.5)*4;
s = [zeros(n-floor((n-m)*(1-degenDegree)),1); rand(floor((n-m)*(1-degenDegree)),1)];

% generate b and c
b = A*x;
c = A'*y+s;

% randperm A and c
perm = randperm(n);
A = A(:,perm);
c = c(perm);

x = x(perm); s = s(perm);

if debug
    fprintf('[m n] = [%d %d]\n',m,n);
    
    if norm(A*x-b)+norm(A'*y+s-c)+sum(x.*s~=0)+sum(x<0)+sum(s<0) < 1e-9
        fprintf('Optimal\n');
        
        if sum(x>0) < m
            fprintf('Primal Degen\n')
        end
        
        if sum(s>0) < n-m
            fprintf('Dual Degen\n');
        end
    else
        fprintf('Not optimal\n');
    end
end
