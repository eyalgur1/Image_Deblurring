function P=LtransExt(X)

[m,n]=size(X);

P{1}=(X(1:m-1,:)-X(2:m,:));
P{2}=(X(:,1:n-1)-X(:,2:n));
P{3}=0.001*X;
c=0.0001;
P{1}=c*P{1};
P{2}=c*P{2};
