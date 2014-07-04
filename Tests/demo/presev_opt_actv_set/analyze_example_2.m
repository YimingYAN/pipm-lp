% VERIFY2
% Original Primal:
% min 2x1+x2
% s.t. x1+x2-x3 = 3
%      x1+2x2-x4 = 4
%      xi>=0 
%      
% which is equivalent to
% min 2x1+x2
% s.t. x1+x2 >= 3
%      x1+2x2 >= 4
%      xi>=0 
%      
% Perturbed Primal
% min (2+lambda)x1 + (1+lambda)x2 + lambda*x3 + lambda*x4
% s.t. x1+x2-x3 = 3
%      x1+2x2-x4 = 4
%      xi >= -lambda
%      
% It's easy to have x3+x4 = 2x1+3x2-7. So perturbed primal is equivalent to
% min (2+3*lambda)x1 + (1+4*lambda)x2 --- slope -(2+3*lambda)/(1+4*lambda)
% s.t. x1+x2 >= 3-lambda  --- slope -1
%      x1+2x2 >= 4-lambda --- slope -1/2
%      x1 >= -lambda      --- slope 0
%      x2 >= -lambda      --- slope infinity
%      
% lambda >= 0
% 
% Potential turnig points: 
% Let -(2+3*lambda)/(1+4*lambda) = infinity >>> no sol
%     -(2+3*lambda)/(1+4*lambda) = 0 >>> no sol
%     -(2+3*lambda)/(1+4*lambda) = -1 >>> lambda = 1
%     -(2+3*lambda)/(1+4*lambda) = -1/2 >>> no sol.
%     
% So the only turning point is lambda=1

close all;
clear;

load verify2.mat;


%%%%%%%%%%%%%%%%%%%%%%
lambda = 1;
%%%%%%%%%%%%%%%%%%%%%%


subplot(1,2,1);
% plot the feasible region for the original problem
hold on;
x1 = 0:0.1:4;

x2 = 3 - x1;
plot(x1,x2);

x2 = 2 - 0.5*x1;
plot(x1,x2);

% objective
% for i=0:0.5:2
%     plot(x1,x1*(-4)/3 + i*lambda,'color','k');
% end
plot(x1,-2*x1,'color','k');
% bounds
x2 = 0 + 0*x1;
plot(x1,x2,'color', 'r', 'linewidth',3);

x1 = 0 + 0*(0:0.1:8);
plot(x1, 0:0.1:8,'color', 'r', 'linewidth',3);

grid on
hold off

title('Original P');
xlabel('x_1');
ylabel('x_2');

subplot(1,2,2);
% plot the feasible region for the perturbed problem
hold on;
x1 = -lambda:0.1:4+lambda;

x2 = 3 - lambda - x1;
plot(x1,x2);

x2 = (4-lambda)/2 - 0.5*x1;
plot(x1,x2);

% objective
% for i=0:1:10
%         plot(x1,x1*(4+lambda)/(4*lambda-3) + i,'color','k');
% end
plot(x1,-x1*(2+3*lambda)/(1+4*lambda),'color','k');
% boudns
x2 = 0 + 0*x1;
plot(x1,x2,'--r', 'linewidth',1);

x1 = 0 + 0*(-lambda:0.1:8);
plot(x1, -lambda:0.1:8,'--r', 'linewidth',1);

x1 = -lambda:0.1:4+lambda;

x2 = -lambda + 0*x1;
plot(x1,x2,'color', 'r', 'linewidth',3);

x1 = -lambda + 0*(-lambda:0.1:8);
plot(x1, -lambda:0.1:8,'color', 'r', 'linewidth',3);

grid on
hold off

title('Perturbed P');
xlabel('x_1');
ylabel('x_2');