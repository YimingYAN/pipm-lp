close all;
clear

n=4;

location = 0;
for itr=1:n
    
    lambda = 1e-2*10^itr;
    
    % primal variables x1 x2
    subplot(n,2,itr+location);
    hold on
    x = -lambda:0.1:1+2*lambda;
    t = -lambda+0*x;
    plot(x,t,'color','red','LineWidth',1)
    plot(t,x,'color','red','LineWidth',1)
    x1=-lambda:0.01:1+lambda;
    x2=(1-x1)/2;
    plot(x1,x2);
    plot(x1,-x1,'color','g');
    for i=-1:0.5:1
        plot(x1,-x1+lambda*i,'color','k');
    end
    hold off
    
    text = ['Primal variables, \lambda = ' num2str(lambda)];
    title(text);
    xlabel('x_1');
    ylabel('x_2');
    
    % dual variables s1 s2
    subplot(n,2,itr+1+location);
    hold on
    s = -lambda:0.1:10;
    t = -lambda+0*s;
    plot(s,t,'color','red','LineWidth',1);
    plot(t,s,'color','red','LineWidth',1);
    s1=(1-lambda)/2:0.01:10;
    s2=-1+2*s1;
    plot(s1,s2);
    plot(s1,min(s2),'color','g');
    grid on
    hold off
    
    text = ['Dual variables, \lambda = ' num2str(lambda)];
    title(text);
    xlabel('s_1');
    ylabel('s_2');
    
    location = location+1;
end

