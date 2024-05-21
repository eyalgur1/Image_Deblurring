%% Hyper-Parameters

R=50;
N=30000;
F_iter_s=[1,10,100,1000];
F_iter_r=10;
PG_iter_s=[10,100,1000];

% LoadDir of output data

LoadDirOut='C:\Users\eyal.gur.STAFF\OneDrive - Technion\06 - Papers\03_NSNAM__Accelerated_Nested_Algorithms_for_Non_Convex_and_Non_Smooth_Optimization_Problems\output\R'; % lab
%LoadDirOut='C:\Users\Gur Eyal\OneDrive - Technion\06 - NSNAM Paper\output\R'; % laptop
%LoadDirOut='D:\Documents\OneDrive - Technion\06 - NSNAM Paper\output\R'; % home

% LoadDir of input data
LoadDirIn='C:\Users\eyal.gur.STAFF\OneDrive - Technion\06 - Papers\03_NSNAM__Accelerated_Nested_Algorithms_for_Non_Convex_and_Non_Smooth_Optimization_Problems\input\R'; % lab
%LoadDirIn='C:\Users\Gur Eyal\OneDrive - Technion\06 - NSNAM Paper\input\R'; % laptop
%LoadDirIn='D:\Documents\OneDrive - Technion\06 - NSNAM Paper\input\R'; % home

load([LoadDirIn,'\hyper_parameters.mat']); % load hyper-parameters

col={'-o', '-s', '-^', '-v', '-d', '-s', '-^', '-v', '-d'}; % markers for plotting
ticks_to_plot=100; % number of jump for marker plotting
legend_plot={'SPA','NAM-F, $s=1$','NAM-F, $s=10$','NAM-F, $s=100$','NAM-F, $s=1000$','ECR-PG, $s=10$','ECR-PG, $s=100$','ECR-PG, $s=1000$'}; % methods
%legend_plot={'SPA','NAM-AG, $s=1$','NAM-AG, $s=10$','NAM-AG, $s=100$','NAM-AG, $s=1000$'}; % methods

%% funv real calculations (not required anymore)
% % setting pars (according to the specification of the tests) for funv real calculation at Bobs and Xor
%
% % reg_g=@(X,alpha)norm(X,1); % l1
% % reg_g=@(X,alpha)norm(X,2)^2; % squared l2
%   reg_g=@(X,alpha)(1-alpha)*norm(X,2)^2+alpha*norm(X,1); % elastic net for alpha\in[0,1]
%
% reg_L=@(X)eye(size(X)); funv_real=zeros(R,1);
% clear pars;
% pars.regfun=reg_g; pars.regtran=reg_L; pars.lambda=1; pars.BC='periodic'; pars.alpha=0.5;
%
% for j=1:R
%     D=([LoadDirIn,'\',num2str(j)]);
%     load([D,'\Bobs.mat']);
%     funv_real(j)=obj_value_real(hyper_parameters.P,hyper_parameters.Xor,Bobs,hyper_parameters.u_real,hyper_parameters.sigma_w,hyper_parameters.sigma_e,pars,hyper_parameters.center);
% end
% funv_real_avg=sum(funv_real)/R;
% funv_real_noNoise=obj_value_real(hyper_parameters.P,hyper_parameters.Xor,hyper_parameters.B,hyper_parameters.u_real,hyper_parameters.sigma_w,hyper_parameters.sigma_e,pars,hyper_parameters.center);

%% Function Values

funv=zeros(N+1,1+length(F_iter_s)+length(PG_iter_s));

for t=1:R
    LoadDir_t=([LoadDirOut,'\',num2str(t)]);
    
    SPA=load([LoadDir_t,'\SPA.mat']);
    funv(:,1)=funv(:,1)+SPA.SPA.funv/R;
    
    for l=1:length(F_iter_s)
        s=F_iter_s(l);
        r=F_iter_r;
        NAMAG=load([LoadDir_t,'\NAMAG_s',num2str(s),'_r',num2str(r),'.mat']);
        funv(:,l+1)=funv(:,l+1)+NAMAG.NAMAG.funv(1:N+1)/R;
    end
    
    for l=1:length(PG_iter_s)
        s=PG_iter_s(l);
        prox=load([LoadDir_t,'\prox_s',num2str(s),'.mat']);
        funv(:,l+1+length(F_iter_s))=funv(:,l+1+length(F_iter_s))+prox.prox.funv(1:N+1)/R;
    end
end

%% Plot Function Values (real)

iter_to_plot=15000;
fh=figure(1);
m=min(min(funv));
tick_to_graph=15;

xplot=1:tick_to_graph:iter_to_plot+1;
for i=1:1+length(F_iter_s)+length(PG_iter_s)
    h=semilogy(funv(xplot,i)-m,col{i}','MarkerIndices',1:ticks_to_plot:length(xplot),'LineWidth',2.5);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    if i>=6
        set(h, 'LineStyle', ':');
    end
    set(gca,'FontSize',20)
    hold on
end
%semilogy(funv_real_noNoise*ones(iter_to_plot+1,1),'black','LineWidth',1.5)
hold off

grid on
xlim([0 length(xplot)])
set(gca,'XTickLabel',strsplit(num2str(xticks)))
xlabel('$\mathbf{Total\ Iterations}$','FontSize',24,'Interpreter','latex')
ylabel('$\mathbf{Function\ Value}\ \left(F_n\right)$','FontSize',24,'Interpreter','latex')
%title('$\mathbf{Average\ Function\ Values\ Over\ 100\ Monte\ Carlo\ Trials}$','FontSize',14,'Interpreter','latex')
legend(legend_plot,'FontSize',24,'Interpreter','latex')
hold off
fh.WindowState = 'maximized';
ax = gca;
exportgraphics(ax,'funv_real.png')

%% Plot Function Values (convegence)

iter_to_plot=15000;
fh=figure(2);
tick_to_graph=15;

xplot=1:tick_to_graph:iter_to_plot+1;
for i=1:1+length(F_iter_s)+length(PG_iter_s)
    h=semilogy(funv(xplot,i)-funv(end,i),col{i},'MarkerIndices',1:ticks_to_plot:length(xplot),'LineWidth',2.5);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    if i>=6
        set(h, 'LineStyle', ':');
    end
    set(gca,'FontSize',20)
    hold on
end
hold off

grid on
xlim([0 length(xplot)])
set(gca,'XTickLabel',strsplit(num2str(xticks)))
xlabel('$\mathbf{Total\ Iterations}$','FontSize',24,'Interpreter','latex')
ylabel('$\Delta F_n$','FontSize',24,'Interpreter','latex')
%title('$\mathbf{Average\ Distance\ of\ Function\ Values\ From\ }F\left(\mathbf{z}^{\mathrm{con}},\mathbf{u}^N\right)\mathbf{\ Over\ 100\ Monte\ Carlo\ Trials}$','FontSize',14,'Interpreter','latex')
legend(legend_plot,'FontSize',24,'Interpreter','latex')
hold off
fh.WindowState = 'maximized';
ax = gca;
exportgraphics(ax,'deltaf.png')

%% Relative Errors (real)

err_real=zeros(N+1,1+length(F_iter_s)+length(PG_iter_s));

for t=1:R
    t
    LoadDir_t=([LoadDirOut,'\',num2str(t)]);
    
    SPA=load([LoadDir_t,'\SPA.mat']);
    err_real(1:N+1,1)=err_real(1:N+1,1)+SPA.SPA.err(1:N+1)/R;
    
    for l=1:length(F_iter_s)
        s=F_iter_s(l);
        r=F_iter_r;
        NAMAG=load([LoadDir_t,'\NAMAG_s',num2str(s),'_r',num2str(r),'.mat']);
        err_real(1:N+1,l+1)=err_real(1:N+1,l+1)+NAMAG.NAMAG.err(1:N+1)/R;
    end
    
    for l=1:length(PG_iter_s)
        s=PG_iter_s(l);
        prox=load([LoadDir_t,'\prox_s',num2str(s),'.mat']);
        err_real(:,l+1+length(F_iter_s))=err_real(:,l+1+length(F_iter_s))+prox.prox.err(1:N+1)/R;
    end
end

%% Plot Relative Errors (real)

iter_to_plot=15000;
figure(3)
m=min(min(err_real));

xplot=1:N+1;
for i=1:1+length(F_iter_s)+length(PG_iter_s)
    h=semilogy(err_real(:,i),col{i}','MarkerIndices',xplot(1:iter_to_plot/ticks_to_plot:end),'LineWidth',2.5);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    if i>=6
        set(h, 'LineStyle', ':');
    end
    set(gca,'FontSize',20)
    hold on
end
hold off

grid on
xlim([0 iter_to_plot+1])
set(gca,'XTickLabel',strsplit(num2str(xticks)))
xlabel('$\mathbf{Iterations}$','FontSize',24,'Interpreter','latex')
ylabel('$\mathbf{Relative\ Error\ (z^{real})}$','FontSize',24,'Interpreter','latex')
%title('$\mathbf{Average\ Relative\ Error\ Over\ 100\ Monte\ Carlo\ Trials\ with\ Respect\ to\ z^{real}}$','FontSize',14,'Interpreter','latex')
legend(legend_plot,'FontSize',24,'Interpreter','latex')
hold off

%% Relative Errors (convergence)

iter_to_plot=15000;
err_con=zeros(iter_to_plot+1,1+length(F_iter_s)+length(PG_iter_s));

err_con_method=zeros(iter_to_plot+1,1);
for t=1:R
    t
    LoadDir_t=([LoadDirOut,'\',num2str(t)]);
    
    SPA=load([LoadDir_t,'\SPA_con.mat']);
    err_con_method=err_con_method+SPA.SPA_con.err/R;
end
err_con(:,1)=err_con_method;

for l=1:length(F_iter_s)
    l
    err_con_method=zeros(iter_to_plot+1,1);
    for t=1:R
        s=F_iter_s(l);
        r=F_iter_r;
        NAMAG=load([LoadDir_t,'\NAMAG_con_s',num2str(s),'_r',num2str(r),'.mat']);
        err_con_method=err_con_method+NAMAG.NAMAG_con.err(1:iter_to_plot+1)/R;
    end
    err_con(:,l+1)=err_con_method;
end

for l=1:length(PG_iter_s)
    l
    err_con_method=zeros(iter_to_plot+1,1);
    for t=1:R
        s=PG_iter_s(l);
        prox=load([LoadDir_t,'\prox_con_s',num2str(s),'.mat']);
        err_con_method=err_con_method+prox.prox_con.err(1:iter_to_plot+1)/R;
    end
    err_con(:,l+length(F_iter_s)+1)=err_con_method;
end

%% Plot Relative Errors (convergence)

iter_to_plot=15000;
fh=figure(4);

xplot=1:iter_to_plot+1;
for i=1:1+length(F_iter_s)+length(PG_iter_s)
    h=semilogy(err_con(:,i),col{i}','MarkerIndices',xplot(1:iter_to_plot/ticks_to_plot:end),'LineWidth',2.5);
    set(h,'MarkerFaceColor',get(h,'Color'));
    if i>=6
        set(h, 'LineStyle', ':');
    end
    set(gca,'FontSize',20)
    hold on
end
hold off

grid on
xlim([0 iter_to_plot+1])
set(gca,'XTickLabel',strsplit(num2str(xticks)))
xlabel('$\mathbf{Total\ Iterations}$','FontSize',24,'Interpreter','latex')
ylabel('$\Delta\mathbf{z}^n$','FontSize',24,'Interpreter','latex')
%title('$\mathbf{Average\ Relative\ Error\ Over\ 100\ Monte\ Carlo\ Trials\ with\ Respect\ to\ z^{con}}$','FontSize',24,'Interpreter','latex')
legend(legend_plot,'FontSize',24,'Interpreter','latex')
hold off
fh.WindowState = 'maximized';
ax = gca;
exportgraphics(ax,'err.png')

%% Plot and Save Real and Blurred Image of Trial #1

X=double(imread('gantrycrane.png'));
plot_image(X(1:256,:,:)/255,5,'$\mathbf{Real\ Image}$')
ax = gca;
exportgraphics(ax,'real_image.png')

load([LoadDirIn,'\1\Bobs.mat']);
plot_image(Bobs(1:256,:,:),6,'$\mathbf{Blurred\ and\ Noisy\ Image}$')
ax = gca;
exportgraphics(ax,'blur_image.png')

%% Plot Output Images of Trial #1

% CRETAE THE IMAGE FILE FIRST

% images_to_plot=[1,2,11]; % number of iterations to plot
% load([LoadDirOut,'\1\SPA_con.mat']);
% 
 fh=figure(7);

% mon=cell(length(images_to_plot),5); % columns - methods, rows - iterations
% for j=images_to_plot
%     ind=find(images_to_plot==j);
%     mon{ind,1}=SPA_con.x(:,:,:,j);
%     load([LoadDirOut,'\1\NAMAG_con_s1_r10.mat']);
%     mon{ind,2}=NAMAG_con.x(:,:,:,j);
%     load([LoadDirOut,'\1\NAMAG_con_s10_r10.mat']);
%     mon{ind,3}=NAMAG_con.x(:,:,:,j);
%     load([LoadDirOut,'\1\NAMAG_con_s100_r10.mat']);
%     mon{ind,4}=NAMAG_con.x(:,:,:,j);
%     load([LoadDirOut,'\1\NAMAG_con_s1000_r10.mat']);
%     mon{ind,5}=NAMAG_con.x(:,:,:,j);
%     
%     load([LoadDirOut,'\1\prox_con_s10.mat']);
%     mon{ind,6}=prox_con.x(:,:,:,j);
%     load([LoadDirOut,'\1\prox_con_s100.mat']);
%     mon{ind,7}=prox_con.x(:,:,:,j);
%     load([LoadDirOut,'\1\prox_con_s1000.mat']);
%     mon{ind,8}=prox_con.x(:,:,:,j);
% end

montage(mon','Size',[3 4],'BorderSize',[2 2],'BackgroundColor','w') % columns - methods, rows - iterations
%pbaspect([5 5 1])
fh.WindowState = 'maximized';
%ax = gca;
%exportgraphics(ax,'output_images.eps')
