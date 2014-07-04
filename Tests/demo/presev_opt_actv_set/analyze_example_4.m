function analyze_example_4
% VERIFY4
% Original P:
% min 4x1+3x2
% s.t. 2x1 +  x2 + x3           = 8
%     -3x1 + 2x2      + x4      = 6
%       x1 + 2x2           + x5 = 6
%
% xi>=0
%
% which is equalvant to
%
% min 4x1+3x2
% s.t. 2x1 +  x2 <= 8
%     -3x1 + 2x2 <= 6
%       x1 + 2x2 <= 6
%
% xi>=0
%
% Perturbed P:
% min (4+lambda)x1+(3+lambda)x2+lambda*x3+lambda*x4+lambda*x5
% s.t. 2x1 +  x2 + x3           = 8
%     -3x1 + 2x2      + x4      = 6
%       x1 + 2x2           + x5 = 6
%
% xi >= -lambda
%
% which is equalvant to
%
% min (4+lambda)x1+(3-4*lambda)x2+20*lambda    
% s.t. 2x1 +  x2 <= 8+lambda                           
%     -3x1 + 2x2 <= 6+lambda                               
%       x1 + 2x2 <= 6+lambda                                 
%
% xi>= -lambda                                               
% note: x3+x4+x5 = 20-5x2;
%
% <=>     
% min (4+lambda)x1+(3-4*lambda)x2 --- slope (4+lambda)/(4*lambda-3)
% s.t. 2x1 +  x2 <= 8+lambda --- slope -2    
%      -3x1 + 2x2 <= 6+lambda --- slope 3/2
%      x1 + 2x2 <= 6+lambda --- slope -1/2
%      x1>= -lambda --- x2 = -lambda, slope 0
%      x2>= -lambda --- x1 = -lambda slope infinity
%
% lambda >= 0
% (4+lambda)/(4*lambda-3) = -2 >>> 2/9
% (4+lambda)/(4*lambda-3) = 3/2 >>> 17/10
% (4+lambda)/(4*lambda-3) = -1/2 >>> no sol
% (4+lambda)/(4*lambda-3) = 0 >>> no sol
% (4+lambda)/(4*lambda-3) = infinity >>> 3/4
%
% So possible turning points are 2/9 3/4 17/10.
% 
% When lambda < 2/9, the optimal solution is (-lambda,-lambda), 
% 2x1 +  x2 = 8+lambda is not an active constraint. So even when lambda =
% 2/9, the optimal solution set won't change.
%
% When lambda = 3/4, the slope of the objective function euqals the slope
% of the constarint x1 = -lambda, which is the active set of the current
% optimal solution (-lambda,-lambda).
%
% When 17/10 > lambda > 3/4, the optimal solution is (-lambda,3-lambda),
% and the active constraints are x1 >= -lambda and -3x1 + 2x2 <= 6+lambda.
% So when lambda = 17/10, the slope of the objective function is the same
% as one of the active constraint, the problem becomes nodegenerate again 
% (turns again).
%
%
% Turning points: all values of lambda that make the slope of 
% the objective function of perturbed problem equal the slope of one of 
% the active constraint.  
%
% When we increase lambda from 0+, the slope of the objective keep 
% increasing from negative to 0. The turning point appears when 
% lambda = 0.75 which makes 4*lambda-3=0, i.e. the slope is infinite,
% which is the same as x1=0. When lambda=0.75, the perturbed primal
% problem has multiple solutions (x2 \in (-0.75,?)), the perturbed dual
% problem becomes degenerate. If we keep increasing lambda, the slope 
% becomes positive, and the optimal solution change from the southwest 
% corner to northwest corner.

close all;
clear;

% options for axis
XMIN = -2;
XMAX = 6;
YMIN = -2;
YMAX = 8;

load verify4.mat;

% problem data --- general form
A(:,3:end) = [];

% plot the orignal problem
subplot(2,2,1);
optx = [];
plotPerturbProb(A,b,0,XMIN, XMAX, YMIN, YMAX,optx);


% lambda = 3/4
subplot(2,2,2);
optx = [];
plotPerturbProb(A,b,3/4,XMIN, XMAX, YMIN, YMAX,optx);

% lambda = 1
subplot(2,2,3);
optx = [];
plotPerturbProb(A,b,1,XMIN, XMAX, YMIN, YMAX,optx);

% lambda = 17/10
subplot(2,2,4);
optx = [];
plotPerturbProb(A,b,17/10,XMIN, XMAX, YMIN, YMAX,optx);

print -dpng 'verify4'
end
% 
% 
% 
function plotPerturbProb(A,b,lambda,XMIN, XMAX, YMIN, YMAX, optx)
% % plot the feasible region for the perturbed problem
hold on;
plot2DRegion(-A,-b-lambda*ones(3,1),-lambda*ones(2,1));
x1 = -lambda:0.1:4+lambda;

x2 = 8+lambda - 2*x1;
plot(x1,x2);

x2 = (6+lambda)/2 + 1.5*x1;
plot(x1,x2);

x2 = (6+lambda)/2 - 0.5*x1;
plot(x1,x2);

% objective
if 4*lambda-3~=0
    t1 = XMIN:0.1:XMAX;
    plot(t1,t1*(4+lambda)/(4*lambda-3),'color','k','linewidth',1.5);
else
    t2 = YMIN:0.1:YMAX;
    t1 = 0*t2-lambda/2;
    plot(t1,t2,'color','k','linewidth',1.5);
end

% bounds
x2 = 0 + 0*x1;
plot(x1,x2,'--r', 'linewidth',1);

x1 = 0 + 0*(-lambda:0.1:8);
plot(x1, -lambda:0.1:8,'--r', 'linewidth',1);

x1 = XMIN:0.1:XMAX;
x2 = -lambda+0*x1;
plot(x1,x2,'color', 'r', 'linewidth',2);

x2 = YMIN:0.1:YMAX;
x1 = -lambda+0*x2;
plot(x1, x2,'color', 'r', 'linewidth',2);

% plot optx
if nargin == 8 && ~isempty(optx)
    plot(optx(:,1),optx(:,2),'marker','*','markersize',2,'color','b');
end

hold off
grid on
if lambda > 0
    title(['Perturbed Prob, lambda = ' num2str(lambda)]);
else
    title('Original Prob');
end
xlabel('x_1');
ylabel('x_2');
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
