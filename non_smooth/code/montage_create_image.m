%% INITIALIZATION

%====================================
% Setting the Regularizer:
% G := g \circ L where g is a function and L is linear transformation
%---------------------
% setting g:
% reg_g=@(X,alpha)norm(X,1); % l1
% reg_g=@(X,alpha)norm(X,2)^2; % squared l2
reg_g=@(X,alpha)(1-alpha)*norm(X,2)^2+alpha*norm(X,1); % elastic net for alpha\in[0,1]
%---------------------
% setting L:
% reg_L=@(X)Ltrans(X);
% reg_L=@(X)LtransExt(X);
% reg_L=@(X)Wtrans(X);
reg_L=@(X)eye(size(X));
%====================================

%====================================
% Setting the Parameters
clear pars;
%---------------------
% setting the prox function:
% pars.prox=@prox_l1;
% pars.prox=@tv_l1;
% pars.prox=@tv_l1plusl2;
% pars.prox=@wave_l1;
% pars.prox=@(x)x/3;
%---------------------
% setting other parameters:
pars.regfun=reg_g;
pars.regtran=reg_L;
pars.lambda=1;
pars.alpha=0.5; % for elastic net in [0,1]
pars.BC='periodic';
% pars.BC='reflexive';
%====================================

%====================================
% Other Settings:
AG_iter_s=[1,1000]; % for AG inner iterations
PG_iter_s=[1000]; % for AG inner iterations
AG_iter_r=10; % for AG inner iterations
R=50; % number of Monte Carlo trials
%====================================

%% ERROR TO CONVERGENCE POINTS

iterations=[1,500,1500]; % setting the number of total iterations
LoadDirIn='C:\Users\eyal.gur\OneDrive - Technion\06 - NSNAM Paper\input\R'; % lab
SaveDirOut='C:\Users\eyal.gur\OneDrive - Technion\06 - NSNAM Paper\output\R';
%LoadDirIn='C:\Users\Gur Eyal\OneDrive - Technion\06 - NSNAM Paper\input\R'; % laptop
%SaveDirOut='C:\Users\Gur Eyal\OneDrive - Technion\06 - NSNAM Paper\output\R';
%LoadDirIn='D:\Documents\OneDrive - Technion\06 - NSNAM Paper\input\R'; % home
%SaveDirOut='D:\Documents\OneDrive - Technion\06 - NSNAM Paper\output\R';

load([LoadDirIn,'\hyper_parameters.mat']); % load the noise used in create_monte_carlo.m

mon=cell(length(iterations),1+length(AG_iter_s)+length(PG_iter_s));

for i=1:length(iterations)
    N=iterations(i);
    
    for j=1:1
        % load the starting point matrices Pobs and Bobs
        SaveDir=([SaveDirOut,'\',num2str(j)]);
        LoadDir=([LoadDirIn,'\',num2str(j)]);
        Pobs=load([LoadDir,'\Pobs.mat']);
        Bobs=load([LoadDir,'\Bobs.mat']);
        
        % run SPA
        Xcon=load([SaveDir,'\SPA.mat']);
        fprintf('SPA: realization %2d/%2d\n', j, R)
        [x_SPA,y_SPA,funv_SPA,err_SPA]=alg_SPA(Pobs.Pobs,hyper_parameters.Xor,Xcon.SPA.x(:,:,:,end),Bobs.Bobs,...
            hyper_parameters.center,hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars);
        
        mon{i,1}=x_SPA(1:256,:,:,end);
        
        % run NAM-AG
        for l=1:length(AG_iter_s)
            s=AG_iter_s(l);
            r=AG_iter_r;
            Xcon=load([SaveDir,'\NAMAG_s',num2str(s),'_r',num2str(r),'.mat']);
            
            fprintf('NAM-AG, s=%2d, r=%2d: realization %2d/%2d\n', s, r, j, R)
            [x_NAMAG,y_NAMAG,funv_NAMAG,err_NAMAG]=alg_NAM_F(Pobs.Pobs,...
                hyper_parameters.Xor,Xcon.NAMAG.x(:,:,:,end),Bobs.Bobs,hyper_parameters.center,...
                hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars,s,r);
            
            mon{i,l+1}=x_NAMAG(1:256,:,:,end);
        end
        
        % run prox
        for l=1:length(PG_iter_s)
            s=PG_iter_s(l);
            Xcon=load([SaveDir,'\prox_s',num2str(s),'.mat']);
            
            fprintf('prox-grad, s=%2d: realization %2d/%2d\n', s, j, R)
            [x_prox,y_prox,funv_prox,err_prox]=alg_ECR_PG(Pobs.Pobs,...
                hyper_parameters.Xor,Xcon.prox.x(:,:,:,end),Bobs.Bobs,hyper_parameters.center,...
                hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars,s);
            
            mon{i,l+1+length(AG_iter_s)}=x_prox(1:256,:,:,end);
        end
    end
end