function analyze_example_4_turningPoint
% VERIFY4
% Original P:
% min -5x1-4x2
% s.t. 2x1 +  x2 + x3           = 8
%     -3x1 + 2x2      + x4      = 6
%       x1 + 2x2           + x5 = 6
%
% xi>=0
%
% which is equalvant to
%
% min -5x1-4x2
% s.t. 2x1 +  x2 <= 8
%     -3x1 + 2x2 <= 6
%       x1 + 2x2 <= 6
%
% xi>=0
%
% Perturbed P:
% min (-5+lambda)x1+(-4+lambda)x2+lambda*x3+lambda*x-5+lambda*x5
% s.t. 2x1 +  x2 + x3           = 8
%     -3x1 + 2x2      + x4      = 6
%       x1 + 2x2           + x5 = 6
%
% xi >= -lambda
%
% which is equalvant to
%
% min (-5+lambda)x1+(-4-4*lambda)x2+20*lambda    
% s.t. 2x1 +  x2 <= 8+lambda                           
%     -3x1 + 2x2 <= 6+lambda                               
%       x1 + 2x2 <= 6+lambda                                 
%
% xi>= -lambda                                               
% note: x3+x4+x5 = 20-5x2;
%
% <=>     
% min (-5+lambda)x1+(-4-4*lambda)x2 --- slope (-5+lambda)/(4*lambda+4)
% s.t. 2x1 +  x2 <= 8+lambda --- slope -2    
%      -3x1 + 2x2 <= 6+lambda --- slope 3/2
%      x1 + 2x2 <= 6+lambda --- slope -1/2
%      x1>= -lambda --- x2 = -lambda, slope 0
%      x2>= -lambda --- x1 = -lambda slope infinity
%
% lambda >= 0
% (-5+lambda)/(4*lambda+4) = -2 >>> no sol
% (-5+lambda)/(4*lambda+4) = 3/2 >>> no sol
% (-5+lambda)/(4*lambda+4) = -1/2 >>> 1
% (-5+lambda)/(4*lambda+4) = 0 >>> 5
% (-5+lambda)/(4*lambda+4) = infinity >>> no sol
%
% So possible turning points are 1 and 5.
% 
% Turning points: all values of lambda that make the slope of 
% the objective function of perturbed problem equal the slope of one of 
% the active constraint.  
%

close all;
clear;

% options for axis
XMIN = -1.1;
XMAX = 5;
YMIN = -1.1;
YMAX = 4;

A = [2 1 1 0 0; -3 2 0 1 0; 1 2 0 0 1];
b = [8 6 6]';
c = [-5 -4 0 0 0]';

% % solve (A,b,c)
% n  = size(A, 2);   
% lb = zeros(n, 1);  ub = inf*ones(n, 1);
% options = optimset('LargeScale', 'off', 'Algorithm','simplex','Display','off');
% [optx,~,~, ~,dual] = linprog(c,[],[],A,b,lb,ub,[],options);
% opts = dual.lower;
% 
% Iactv_x = optx>1e-16;
% Iactv_s = opts>1e-16;
% 
% lambda_bar = min( [optx(Iactv_x);opts(Iactv_s)])/max(svd(inv(A(:,Iactv_x))*A(:,Iactv_s)))
% lambda_bar/sqrt(2)

% problem data --- general form
A(:,3:end) = [];
lambda = 0:1/3:1;
for i =1:4
% plot the orignal problem
subplot(2,2,i);
optx = [];
plotPerturbProb(A,b,lambda(i),XMIN, XMAX, YMIN, YMAX,optx);
end

print -dpng 'verify4'
end
% 
% 
% 
function plotPerturbProb(A,b,lambda,XMIN, XMAX, YMIN, YMAX, optx)
% % plot the feasible region for the perturbed problem
hold on;
plot2DRegion(-A,-b-lambda*ones(3,1),-lambda*ones(2,1));

x1 = -lambda:0.1:-5+lambda;

x2 = 8+lambda - 2*x1;
plot(x1,x2);

x2 = (6+lambda)/2 + 1.5*x1;
plot(x1,x2);

x2 = (6+lambda)/2 - 0.5*x1;
plot(x1,x2);

% objective
t1 = XMIN:0.1:XMAX;

for const = 0:1:5
plot(t1,t1*(-5+lambda)/(4*lambda+4)+const,'--k','linewidth',0.2);
end


% bounds
x1 = XMIN:0.1:XMAX;
x2 = 0 + 0*x1;
plot(x1,x2,'--r', 'linewidth',1);
x2 = -lambda+0*x1;
plot(x1,x2,'color', 'r', 'linewidth',2);

x2 = YMIN:0.1:YMAX;
x1 = 0 + 0*x2;
plot(x1, x2,'--r', 'linewidth',1);
x1 = -lambda+0*x2;
plot(x1, x2,'color', 'r', 'linewidth',2);


% plot optx
if nargin == 8 && ~isempty(optx)
    plot(optx(:,1),optx(:,2),'marker','*','markersize',2,'color','b');
end

hold off
grid on
if lambda > 0
    title(['Perturbed Prob, \lambda = ' num2str(lambda)]);
else
    title('Original Prob');
end
xlabel('x_{1}');
ylabel('x_{2}');
axis([XMIN XMAX YMIN YMAX])
end

function h = plot2DRegion(A,b,lb)

[m,n]=size(A);
transp=0.5;

if size(b,1)==1
    b=b';
end


if size(lb,1)==1
    lb=lb';
end
A=[A;eye(n)];
b=[b;lb];
m=m+n;


c=[rand rand rand];
% c = 'r';

eq=zeros(2,1);
X=zeros(2,1);
for i=1:(m-1)
    for j=(i+1):m
        try
            x=A([i j],:)\b([i j]);
            if and(min((A*x-b))>-1e-6,min((A*x-b))<Inf)
                X=[X,x];
                eq=[eq,[i j]'];
            end
        end
    end
end


[rad,col]=size(X);
xm=mean(X(:,2:end),2);
Xdiff=X(:,2:end);

for j=1:(col-1)
    Xdiff(:,j)=Xdiff(:,j)-xm;
    Xdiff(:,j)=Xdiff(:,j)/norm(Xdiff(:,j));
end
costhe=zeros((col-1),1);

for j=1:(col-1)
    costhe(j)=Xdiff(:,1)'*Xdiff(:,j);
end

[cc,ind]=min(abs(costhe));
ref2=Xdiff(:,ind(1))-(Xdiff(:,ind(1))'*Xdiff(:,1))*Xdiff(:,1);
ref2=ref2'/norm(ref2);

for j=1:(col-1)
    if ref2*Xdiff(:,j)<0
        costhe(j)=-2-costhe(j);
    end
end
[sooo,ind3]=sort(costhe);
h = patch(X(1,ind3+1)',X(2,ind3+1)',c);
set(h,'FaceAlpha',transp);
end
