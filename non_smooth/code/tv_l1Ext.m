function [prox_sol,i]=tv_l1Ext(X,lambda,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the proximal mapping of g(x)=TV(x)=\norm{Lx}_{1}
% using FISTA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%
% X ....................... the given vector
% lambda................... the parameter of the prox mapping
% N ....................... number of iterations

% OUTPUT:
%
% prox_sol.................  the solution of the prox computation

%Define the Projection onto the box
%l=0;
l=-Inf;
%u=1;
u=Inf;
if((l==-Inf)&(u==Inf))
    project=@(x)x;
elseif (isfinite(l)&(u==Inf))
    project=@(x)(((l<x).*x)+(l*(x<=l)));
elseif (isfinite(u)&(l==-Inf))
     project=@(x)(((x<u).*x)+((x>=u)*u));
elseif ((isfinite(u)&isfinite(l))&(l<u))
    project=@(x)(((l<x)&(x<u)).*x)+((x>=u)*u)+(l*(x<=l));
else
    error('lower and upper bound l,u should satisfy l<u');
end

[m,n]=size(X);
epsilon=1e-4;
prnt=1;

% intial point
P{1}=zeros(m-1,n);
P{2}=zeros(m,n-1);
R{1}=zeros(m-1,n);
R{2}=zeros(m,n-1);
P{3}=zeros(m,n);
R{3}=zeros(m,n);

tk=1;
tkp1=1;
count=0;
i=0;

D=zeros(m,n);
fun_all=[];
fun_pri_all=[];
fval=inf;

while((i<N)&(count<5))
    fold=fval;
    % updating the iteration counter
    i=i+1;

    % Storing the old value of the current solution
    Dold=D;
    Pold=P;
    tk=tkp1;

    %Computing the gradient of the objective function
    D=project(X-(lambda/2)*LforwardExt(R));
    Q=LtransExt(D);

    % Taking a step towards minus of the gradient
    P{1}=R{1}+1/(4.5*lambda)*Q{1};
    P{2}=R{2}+1/(4.5*lambda)*Q{2};
    P{3}=R{3}+1/(4.5*lambda)*Q{3};
    
    % Peforming the projection step
    P{1}=P{1}./(max(abs(P{1}),1));
    P{2}=P{2}./(max(abs(P{2}),1));
    P{3}=P{3}./(max(abs(P{3}),1));
    
    %Updating R and t
    tkp1=(1+sqrt(1+4*tk^2))/2;
    R{1}=P{1}+(tk-1)/(tkp1)*(P{1}-Pold{1});
    R{2}=P{2}+(tk-1)/tkp1*(P{2}-Pold{2});
    R{3}=P{3}+(tk-1)/tkp1*(P{3}-Pold{3});
    
    fval=norm(project(X-(lambda/2)*LforwardExt(P)),'fro')^2;
    fun_all=[fun_all;fval];

%     A2=Ltrans(D);
%     fval_pri=norm(D-X,'fro')^2+lambda*(sum(sum(abs(A2{1})))+sum(sum(abs(A2{2}))));
%     fun_pri_all=[fun_pri_all;fval_pri];

    re=norm(D-Dold,'fro')/norm(D,'fro');
    if (re<epsilon)
        count=count+1;
    else
        count=0;
    end
    if(prnt)
%         fprintf('iter= %5d value = %10.10f %10.10f\n',i,fval_pri,norm(D-Dold,'fro')/norm(D,'fro'));
%         if (fval>fold)
%             fprintf('  *\n');
%         else
%             fprintf('   \n');
%         end
    end
end
prox_sol=D;

end







