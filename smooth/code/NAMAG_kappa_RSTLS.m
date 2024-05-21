function [x_APB,y_APB,val_all,rlerr_all]=NAMAG_kappa_RSTLS(P,Xor,B,center,sigma_w,sigma_e,N,pars,s,r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function seeks a solution to the regularized structued total least
% squares problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%
% P ....................... the noised PSF
% Xor ..................... the original image
% B ....................... the noised image
% center .................. the center of the PSF
% sigma_w.................. the standard deviation of the righthand side vector
% sigma_e ................. the standard deviation of the structure components
% N ....................... the number of iterations
% pars .................... Parameters structure
%   pars.prox ............. Type of regularizer
%   pars.regfun ........... The regularizer
%   pars.regtran .......... The linear transformation
%   pars.BC ............... Type of boundary conditions: 'reflexive' or 'periodic.
%
% OUTPUT:
%
% x_APB ................. the x solution of the RSTLS problem
% y_APB ................. the y solution of the RSTLS problem
% val_all ............... objective function values at each iteration
% rlerr_all ............. error at each iteration

%% INITIALIZATION

rng('default');

switch pars.BC
    case 'reflexive'
        trans=@(X)dct2(X);
        itrans=@(X)idct2(X);
        % computng the eigenvalues of the blurring matrix
        e1=zeros(size(B));
        e1(1,1)=1;
        eigval=@(X)dct2(dctshift(X,center))./dct2(e1);
    case 'periodic'
        trans=@(X)fft2(X);
        itrans=@(X)ifft2(X);
        % computng the eigenvalues of the blurring matrix
        eigval=@(X)fft2(circshift(X,1-center));
    otherwise
        error('Invalid boundary conditions should be reflexive or periodic');
end

% Detect the structure of the PSF matrix P
[p,S]=detect_structure(P);

% Generate starting points
X0=B;
%X0=rand(size(B));
%y0=rand(p,1);
y0=zeros(p,1);

% Pading the PSF S
Sbig=zeros([size(B),p]);
for i=1:p
    Sbig(:,:,i)=padPSF(S(:,:,i),size(B));
end

% Finding the vetor of eigenvalues of the matrix Sbig
DSbig=zeros([size(B),p]);
for i=1:p
    DSbig(:,:,i)=eigval(Sbig(:,:,i));
end

% Computing the two dimensional transform of Bobs
Btrans=trans(B);

% Computing eigenvalues of the first noised PSF
E0=zeros(size(P));
for i=1:p
    E0=E0+y0(i)*S(:,:,i);
end
C0=P+E0;

% Pading the PSF
C0big=padPSF(C0,size(B));

% Finding the vetor of eigenvalues of the matrix Cbig
DC0=eigval(C0big);

% Pading the PSF
Pbig=padPSF(P,size(B));

% Finding the vetor of eigenvalues of the matrix Cbig
DP=eigval(Pbig);

% Initializion of NAMAG
X=X0;
vecX=reshape(X,size(B,1)*size(B,2),1);
DC=DC0;
y=y0;
val_all=[];
rlerr_all=[];
aa=norm(Xor,'fro');

iter=0; % accumulative iteraion counter (inclusing the first u-update)
out_iter=0; % outer iteraion counter (inclusing the first u-update)

%%
while iter<N
    V=X;
    
    out_iter=out_iter+1;
    in_iter=0; % AG inner itearion counter
    
    L=2*max(max(abs(DC).^2))+(sigma_w^2)*(pars.lambda)*2; % Lipschitz constsnt
    kappa=L/((sigma_w^2)*(pars.lambda)*2); % condition number
    AGiter=ceil(log(0.5)/log(1-(1/sqrt(kappa)))); % set the number of inner AG iterations

    %t=1; % stepsize intialization when kappa is not known
    
    %%%%%%%%%%%%%The x-step%%%%%%%%%%%%%
    while in_iter<AGiter
        in_iter=in_iter+1;
        vecX_prev=vecX;
        %t_prev=t; % revious stepsize updtate when kappa is not known
        
        % setting the gradient
        Vtrans=DC.*trans(V);
        G=reshape(2*real(itrans(conj(DC).*(Vtrans-Btrans))),size(B,1)*size(B,2),1);
        vecV=reshape(V,size(B,1)*size(B,2),1);
        G=G+(sigma_w^2)*(pars.lambda)*2*vecV; % the gradient
        
        % gradient step
        vecX=vecV-(1/L)*G;
        
        % AG update rule when kappa is known
        vecV=vecX+((sqrt(kappa)-1)/(sqrt(kappa)+1))*(vecX-vecX_prev);
        
        % AG update rule when kappa is not known
        %t=(1+sqrt(1+4*(t^2)))/2;
        %vecV=vecX+((t_prev-1)/t)*(vecX-vecX_prev);
        
        V=reshape(vecV,size(B,1),size(B,1));
        X=reshape(vecX,size(B,1),size(B,1));
        AX=real(itrans(DC.*trans(X)));
        
        if in_iter<AGiter
            % Compute the objective value
            val_all=[val_all;obj_value(AX,X,B,y,sigma_w,sigma_e,pars)];
        
            % Compute the error
            rlerr_all=[rlerr_all;norm(X-Xor,'fro')/aa];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%The y-step%%%%%%%%%%%%%
    PX=real(itrans(DP.*trans(X)));
    a=vect(PX-B);
    m=size(a);
    
    % Computting the matrix [E_1x,...,E_px]
    E1=zeros(m(1),p);
    XF=trans(X);
    for i=1:p
        E1(:,i)=vect(real(itrans(DSbig(:,:,i).*XF)));
    end
    E2=(sigma_w/sigma_e)^2*eye(p)+E1'*E1;
    y=-E2\(E1'*a);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E=zeros(size(P));
    for i=1:p
        E=E+y(i)*S(:,:,i);
    end
    C=P+E;
    
    % Pading the PSF
    Cbig=padPSF(C,size(B));
    
    % Finding the vetor of eigenvalues of the matrix Cbig
    DC=eigval(Cbig);
    
    AX=real(itrans(DC.*trans(X)));
    
    % Compute the objective value
    val=obj_value(AX,X,B,y,sigma_w,sigma_e,pars);
    val_all=[val_all;val];
    
    % Compute the error
    rlerr=norm(X-Xor,'fro')/aa;
    rlerr_all=[rlerr_all;rlerr];
    
    fprintf('#%d value=%f   RleErr=%f   L=%10.4f   #inner iterations=%d\n',iter,val,rlerr,L,1);
    
    iter=iter+AGiter; % update the accumulative iteration counter
end

%% OUTPUT

rlerr_all=[norm(X0-Xor,'fro')/aa;rlerr_all]; % add error at starting point
AX0=real(itrans(DC0.*trans(X0)));
val_all=[obj_value(AX0,X0,B,y0,sigma_w,sigma_e,pars);val_all]; % add function value at starting point

% Saving the obtained soultion
x_APB=X;
y_APB=y;

end

%% AUXILIARY FUNCTIONS

% Compute the objective value
function v=obj_value(AX,X,B,y,sigma_w,sigma_e,pars)

V=pars.regtran(X);
[~,m]=size(V);
res_reg=0;
cell_bit=iscell(V);
if cell_bit==0
    v=X(:);
    res_reg=pars.regfun(v);
else
    for i=1:m
        VV=V{i};
        vv=VV(:);
        res_reg=res_reg+pars.regfun(vv);
    end
end
v=(sigma_w/sigma_e)^2*norm(y)^2+norm(AX-B,'fro')^2+(sigma_w^2)*pars.lambda*res_reg;
end

% Transform matrix to vector
function x=vect(X)
x=X(:);
end

% Compute the Lip. constant of the gradient of H with respect to variable X using bactracking
function [L,X]=backt(X,G,DC,B,y,pars)
% Initializion of the parameters
L=0.1;
eta=1.1;

% The gradient step
Z=X-(1/L)*G;
% The prox step
Pr=prox_l1(Y,pars.lambda,L);
% Computing the four terms of the Descent Lemma (2 function values (at X and at Pr),
% Inner Product Term and Norm Term
AX=real(idct2(DC.*dct2(X)));
APr=real(idct2(DC.*dct2(Pr)));
FX=obj_value(AX,X,B,y,0);
FPr=obj_value(APr,Pr,B,y,0);
IPT=sum(sum(G.*(Pr-X)));
NT=(L/2)*norm(Pr-X,'fro')^2;
while (FPr>FX+IPT+NT)
    L=L*eta;
    Z=X-(1/L)*G;
    Pr=prox_l1(Z,pars.lambda,L);
    APr=real(idct2(DC.*dct2(Pr)));
    FPr=obj_value(APr,Pr,B,y,0);
    IPT=sum(sum(G.*(Pr-X)));
    NT=(L/2)*norm(Pr-X,'fro')^2;
end
X=Pr;
end