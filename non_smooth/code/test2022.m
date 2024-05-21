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
AG_iter_s=[4000]; % for AG inner iterations
AG_iter_r=2; % for AG inner iterations
R=50; % number of Monte Carlo trials
%====================================

%% ERROR TO CONVERGENCE POINTS

% N=15000; % setting the number of total iterations
% LoadDirIn='C:\Users\eyal.gur\OneDrive - Technion\06 - NSNAM Paper\input\R'; % lab
% SaveDirOut='C:\Users\eyal.gur\OneDrive - Technion\06 - NSNAM Paper\output\R';
% LoadDirIn='C:\Users\Gur Eyal\OneDrive - Technion\06 - NSNAM Paper\input\R'; % laptop
% SaveDirOut='C:\Users\Gur Eyal\OneDrive - Technion\06 - NSNAM Paper\output\R';
% LoadDirIn='D:\Documents\OneDrive - Technion\06 - NSNAM Paper\input\R'; % home
% SaveDirOut='D:\Documents\OneDrive - Technion\06 - NSNAM Paper\output\R';
%
% load([LoadDirIn,'\hyper_parameters.mat']); % load the noise used in create_monte_carlo.m
%
% for j=41:50
%     load the starting point matrices Pobs and Bobs
%     SaveDir=([SaveDirOut,'\',num2str(j)]);
%     LoadDir=([LoadDirIn,'\',num2str(j)]);
%     Pobs=load([LoadDir,'\Pobs.mat']);
%     Bobs=load([LoadDir,'\Bobs.mat']);

%     % run SPA
%     Xcon=load([SaveDir,'\SPA.mat']);
%     fprintf('SPA: realization %2d/%2d\n', j, R)
%     [x_SPA,y_SPA,funv_SPA,err_SPA]=SPA_RSTLS(Pobs.Pobs,hyper_parameters.Xor,Xcon.SPA.x(:,:,:,end),Bobs.Bobs,...
%         hyper_parameters.center,hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars);
%
%         SPA_con=struct; SPA_con.x=x_SPA; SPA_con.y=y_SPA; SPA_con.funv=funv_SPA; SPA_con.err=err_SPA; SPA_con.N=N;
%         save([SaveDir,'/SPA_con.mat'],'SPA_con')
%
%     % run NAM-AG
%         for l=1:length(AG_iter_s)
%             s=AG_iter_s(l);
%             r=AG_iter_r;
%             Xcon=load([SaveDir,'\NAMAG_s',num2str(s),'_r',num2str(r),'.mat']);
%
%             fprintf('NAM-AG, s=%2d, r=%2d: realization %2d/%2d\n', s, r, j, R)
%             [x_NAMAG,y_NAMAG,funv_NAMAG,err_NAMAG]=NAMAG_RSTLS(Pobs.Pobs,...
%                 hyper_parameters.Xor,Xcon.NAMAG.x(:,:,:,end),Bobs.Bobs,hyper_parameters.center,...
%                 hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars,s,r);
%
%             NAMAG_con=struct; NAMAG_con.x=x_NAMAG; NAMAG_con.y=y_NAMAG; NAMAG_con.funv=funv_NAMAG; NAMAG_con.err=err_NAMAG; NAMAG_con.N=N;
%             save([SaveDir,'/NAMAG_con_','s',num2str(s),'_r',num2str(r),'.mat'],'NAMAG_con')
%         end

% run NAM-AG
%     for l=1:length(AG_iter_s)
%         s=AG_iter_s(l);
%         Xcon=load([SaveDir,'\prox_s',num2str(s),'.mat']);
%
%         fprintf('prox-grad, s=%2d: realization %2d/%2d\n', s, j, R)
%         [x_prox,y_prox,funv_prox,err_prox]=prox_RSTLS(Pobs.Pobs,...
%             hyper_parameters.Xor,Xcon.prox.x(:,:,:,end),Bobs.Bobs,hyper_parameters.center,...
%             hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars,s);
%
%         prox_con=struct; prox_con.x=x_prox; prox_con.y=y_prox; prox_con.funv=funv_prox; prox_con.err=err_prox; prox_con.N=N;
%         save([SaveDir,'/prox_con_','s',num2str(s),'.mat'],'prox_con')
%     end
% end

%% ERROR TO REAL POINT

N=4000; % setting the number of total iterations

% setting the directory with initial matrices (created by create_monte_carlo.m)
LoadDirIn='C:\Users\eyal.gur\OneDrive - Technion\06 - NSNAM Paper\input\R'; % lab
SaveDirOut='C:\Users\eyal.gur\OneDrive - Technion\06 - NSNAM Paper\output\R';
%LoadDirIn='C:\Users\Gur Eyal\OneDrive - Technion\06 - NSNAM Paper\input\R'; % laptop
%SaveDirOut='C:\Users\Gur Eyal\OneDrive - Technion\06 - NSNAM Paper\output\R';
%LoadDirIn='D:\Documents\OneDrive - Technion\06 - NSNAM Paper\input\R'; % home
%SaveDirOut='D:\Documents\OneDrive - Technion\06 - NSNAM Paper\output\R';

load([LoadDirIn,'\hyper_parameters.mat']); % load the noise used in create_monte_carlo.m

for j=1:1
    % load the starting point matrices Pobs and Bobs
    SaveDir=([SaveDirOut,'\',num2str(j)]);
    LoadDir=([LoadDirIn,'\',num2str(j)]);
    Pobs=load([LoadDir,'\Pobs.mat']);
    Bobs=load([LoadDir,'\Bobs.mat']);
    
    % run SPA
    %     fprintf('SPA: realization %2d/%2d\n', j, R)
    %     [x_SPA,y_SPA,funv_SPA,err_SPA]=SPA_RSTLS(Pobs.Pobs,hyper_parameters.Xor,hyper_parameters.Xor,Bobs.Bobs,...
    %         hyper_parameters.center,hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars);
    %
    %     SPA=struct; SPA.x=x_SPA; SPA.y=y_SPA; SPA.funv=funv_SPA; SPA.err=err_SPA; SPA.N=N;
    %save([SaveDir,'/SPA.mat'],'SPA')
    
    % run NAM-F
    for l=1:length(AG_iter_s)
        s=AG_iter_s(l);
        r=AG_iter_r;
        
        fprintf('NAM-FR, s=%2d, r=%2d: realization %2d/%2d\n', s, r, j, R)
        [x_NAMFR,y_NAMFR,funv_NAMFR,err_NAMFR]=alg_NAM_FR(Pobs.Pobs,...
            hyper_parameters.Xor,hyper_parameters.Xor,Bobs.Bobs,hyper_parameters.center,...
            hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars,s,r);
        
        NAMFR=struct; NAMFR.x=x_NAMFR; NAMFR.y=y_NAMFR; NAMFR.funv=funv_NAMFR; NAMFR.err=err_NAMFR; NAMFR.N=N;
        %save([SaveDir,'/NAMFR_','s',num2str(s),'_r',num2str(r),'.mat'],'NAMFR')
    end
    
    % run NAM-F
    for l=1:length(AG_iter_s)
        s=AG_iter_s(l);
        r=AG_iter_r;
        
        fprintf('NAM-F, s=%2d, r=%2d: realization %2d/%2d\n', s, r, j, R)
        [x_NAMF,y_NAMF,funv_NAMF,err_NAMF]=alg_NAM_F(Pobs.Pobs,...
            hyper_parameters.Xor,hyper_parameters.Xor,Bobs.Bobs,hyper_parameters.center,...
            hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars,s,r);
        
        NAMF=struct; NAMF.x=x_NAMF; NAMF.y=y_NAMF; NAMF.funv=funv_NAMF; NAMF.err=err_NAMF; NAMF.N=N;
        %save([SaveDir,'/NAMF_','s',num2str(s),'_r',num2str(r),'.mat'],'NAMF')
    end
    
    % run ECR-PG
    %     for l=1:length(AG_iter_s)
    %         s=AG_iter_s(l);
    %
    %         fprintf('prox-grad, s=%2d: realization %2d/%2d\n', s, j, R)
    %         [x_prox,y_prox,funv_prox,err_prox]=alg_ECR_PG(Pobs.Pobs,...
    %             hyper_parameters.Xor,hyper_parameters.Xor,Bobs.Bobs,hyper_parameters.center,...
    %             hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars,s);
    %
    %         prox=struct; prox.x=x_prox; prox.y=y_prox; prox.funv=funv_prox; prox.err=err_prox; prox.N=N;
    %         %save([SaveDir,'/prox_','s',num2str(s),'.mat'],'prox')
    %     end
end