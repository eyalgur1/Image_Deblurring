function [prox_sol,i]=wave_l1(X,lam,t,N)
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

lambda=lam/t;

[m,n]=size(X);
epsilon=1e-4;
prnt=1;

% intial point
P=zeros(size(X));
R=zeros(size(X));

tk=1;
tkp1=1;
count=0;
i=0;

D=zeros(size(X));
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

%     %Computing the gradient of the objective function
%     D=project(X-(lambda/2)*Wforward(R));
%     Q=Wtrans(D);

    % Taking a step towards minus of the gradient
    [L,P,D]=backt(R,X,lambda);
 
    %Updating R and t
    tkp1=(1+sqrt(1+4*tk^2))/2;
    R=P+(tk-1)/(tkp1)*(P-Pold);

    fval=norm(project(X-(lambda/2)*Wforward(P)),'fro')^2;
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


% Compute the Lip. constant of the gradient of H with respect to variable X using bactracking
function [L,P,D]=backt(R,X,lambda)
    % Initializion of the parameters
    L=0.1;
    eta=1.1;
    
    % The gradient step
    D=X-(lambda/2)*Wforward(R);
    G=-lambda*Wtrans(D);
    Y=R-(1/L)*G;
    % The prox step
    Pr=Y./(max(abs(Y),1));
    % Computing the four terms of the Descent Lemma (2 function values (at X and at Pr),
    % Inner Product Term and Norm Term
    FX=norm(X-(lambda/2)*Wforward(R),'fro')^2;
    FPr=norm(X-(lambda/2)*Wforward(Pr),'fro')^2;
    IPT=sum(sum(G.*(Pr-R)));
    NT=(L/2)*norm(Pr-R,'fro')^2;
    while (FPr>FX+IPT+NT)
        L=L*eta;
        Y=R-(1/L)*G;
        Pr=Y./(max(abs(Y),1));
        FPr=norm(X-(lambda/2)*Wforward(Pr),'fro')^2;
        IPT=sum(sum(G.*(Pr-R)));
        NT=(L/2)*norm(Pr-R,'fro')^2;
    end
    P=Pr;
end






