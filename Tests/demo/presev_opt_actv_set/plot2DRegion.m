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