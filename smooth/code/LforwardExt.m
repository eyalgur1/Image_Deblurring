function X=LforwardExt(P)

[m2,n2]=size(P{1});
[m1,n1]=size(P{2});
[m3,n3]=size(P{3});



if (n2~=n1+1)
    error('dimensions are not consistent')
end
if(m1~=m2+1)
    error('dimensions are not consistent')
end

m=m2+1;
n=n2;
c=0.0001;
P{1}=c*P{1};
P{2}=c*P{2};
X=0.001*P{3};
X(1:m-1,:)=P{1};
X(:,1:n-1)=X(:,1:n-1)+P{2};
X(2:m,:)=X(2:m,:)-P{1};
X(:,2:n)=X(:,2:n)-P{2};


