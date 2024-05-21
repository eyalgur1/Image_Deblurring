%% INITIALIZATION

%====================================
% Setting the Regularizer:
% G := g \circ L where g is a function and L is linear transformation
%---------------------
% setting g:
% reg_g=@(X)norm(X,1);
reg_g=@(X)norm(X,2)^2;
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
pars.lambda=0.02;
pars.BC='periodic';
% pars.BC='reflexive';
%====================================

%====================================
% Other Settings:
AG_iter_s=[1,10,50,100,200]; % for AG inner iterations
AG_iter_r=10; % for AG inner iterations
R=100; % number of Monte Carlo trials
%====================================

%% ERROR TO CONVERGENCE POINTS

N=1; % setting the number of total iterations
%LoadDir='C:\Users\eyal.gur\OneDrive - Technion\04 - NAM Paper\14 - Numerical Analysis - 27.6.21\R'; % setting the directory with initial matrices (created by create_monte_carlo.m)
%LoadDir='C:\Users\Gur Eyal\OneDrive - Technion\04 - NAM Paper\14 - Numerical Analysis - 27.6.21\R';
LoadDir='D:\Documents\OneDrive - Technion\04 - NAM Paper\14 - Numerical Analysis - 27.6.21\R';
load([LoadDir,'\hyper_parameters.mat']); % load the noise used in create_monte_carlo.m

for j=1:R
    % load the starting point matrices Pobs and Bobs
    SaveDir=([LoadDir,'\',num2str(j)]);
    Pobs=load([SaveDir,'\Pobs.mat']);
    Bobs=load([SaveDir,'\Bobs.mat']);
    
    % run SPA
    Xcon=load([SaveDir,'\SPA.mat']);
    fprintf('SPA: realization %2d/%2d\n', j, R)
    [x_SPA,y_SPA,funv_SPA,err_SPA]=SPA_RSTLS(Pobs.Pobs,hyper_parameters.Xor,Xcon.SPA.x,Bobs.Bobs,...
        hyper_parameters.center,hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars);
    
        SPA_result=struct; SPA_result.x=x_SPA; SPA_result.y=y_SPA; SPA_result.funv=funv_SPA; SPA_result.err=err_SPA; SPA_result.N=N;
        %save([SaveDir,'/SPA_result.mat'],'SPA_result')
    
    % run NAM-AG
        for l=1:length(AG_iter_s)
            s=AG_iter_s(l);
            r=AG_iter_r;
            Xcon=load([SaveDir,'\NAMAG_s',num2str(s),'_r',num2str(r),'.mat']);
    
            fprintf('NAM-AG, s=%2d, r=%2d: realization %2d/%2d\n', s, r, j, R)
            [x_NAMAG,y_NAMAG,funv_NAMAG,err_NAMAG]=NAMAG_RSTLS(Pobs.Pobs,...
                hyper_parameters.Xor,Xcon.NAMAG.x,Bobs.Bobs,hyper_parameters.center,...
                hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars,s,r);
    
            NAMAG_result=struct; NAMAG_result.x=x_NAMAG; NAMAG_result.y=y_NAMAG; NAMAG_result.funv=funv_NAMAG; NAMAG_result.err=err_NAMAG; NAMAG_result.N=N;
            %save([SaveDir,'/NAMAG_result_','s',num2str(s),'_r',num2str(r),'.mat'],'NAMAG_result')
        end
end

%% ERROR TO REAL POINT

% N=50000; % setting the number of total iterations
% LoadDir='C:\Users\eyal.gur\OneDrive - Technion\04 - NAM Paper\14 - Numerical Analysis - 27.6.21\R'; % setting the directory with initial matrices (created by create_monte_carlo.m)
% LoadDir='C:\Users\Gur Eyal\OneDrive - Technion\04 - NAM Paper\14 - Numerical Analysis - 27.6.21\R';
% LoadDir='D:\Documents\OneDrive - Technion\04 - NAM Paper\14 - Numerical Analysis - 27.6.21\R';
% load([LoadDir,'\hyper_parameters.mat']); % load the noise used in create_monte_carlo.m
%
% for j=1:R
%     load the starting point matrices Pobs and Bobs
%     SaveDir=([LoadDir,'\',num2str(j)]);
%     Pobs=load([SaveDir,'\Pobs.mat']);
%     Bobs=load([SaveDir,'\Bobs.mat']);
%
%     run SPA
%     fprintf('SPA: realization %2d/%2d\n', j, R)
%     [x_SPA,y_SPA,funv_SPA,err_SPA]=SPA_RSTLS(Pobs.Pobs,hyper_parameters.Xor,hyper_parameters.Xor,Bobs.Bobs,...
%         hyper_parameters.center,hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars);
%
%     SPA=struct; SPA.x=x_SPA; SPA.y=y_SPA; SPA.funv=funv_SPA; SPA.err=err_SPA; SPA.N=N;
%     save([SaveDir,'/SPAt.mat'],'SPA')
%
%     run NAM-AG
%     for l=1:length(AG_iter_s)
%         s=AG_iter_s(l);
%         r=AG_iter_r;
%
%         fprintf('NAM-AG, s=%2d, r=%2d: realization %2d/%2d\n', s, r, j, R)
%         [x_NAMAG,y_NAMAG,funv_NAMAG,err_NAMAG]=NAMAG_RSTLS(Pobs.Pobs,...
%             hyper_parameters.Xor,hyper_parameters.Xor,Bobs.Bobs,hyper_parameters.center,...
%             hyper_parameters.sigma_w,hyper_parameters.sigma_e,N,pars,s,r);
%
%         NAMAG=struct; NAMAG.x=x_NAMAG; NAMAG.y=y_NAMAG; NAMAG.funv=funv_NAMAG; NAMAG.err=err_NAMAG; NAMAG.N=N;
%         save([SaveDir,'/NAMAG_','s',num2str(s),'_r',num2str(r),'.mat'],'NAMAG')
%     end
% end