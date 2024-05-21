%% Hyper-Parameters

R=100;
N=50000;
AG_iter_s=[1,10,50,100,200];
AG_iter_r=10;
%LoadDir='C:\Users\Gur Eyal\OneDrive - Technion\04 - NAM Paper\14 - Numerical Analysis - 27.6.21\R';
LoadDir='C:\Users\eyal.gur\OneDrive - Technion\04 - NAM Paper\14 - Numerical Analysis - 27.6.21\R';
%LoadDir='D:\Documents\OneDrive - Technion\04 - NAM Paper\14 - Numerical Analysis - 27.6.21\R';
load([LoadDir,'\hyper_parameters.mat']);
funv_real=zeros(R,1);
reg_g=@(X)norm(X,2)^2;
reg_L=@(X)eye(size(X));
clear pars;
pars.regfun=reg_g; pars.regtran=reg_L; pars.lambda=0.02; pars.BC='periodic';
for j=1:R
    D=([LoadDir,'\',num2str(j)]);
    load([D,'\Bobs.mat']);
    funv_real(j)=obj_value_real(hyper_parameters.P,hyper_parameters.Xor,Bobs,hyper_parameters.u_real,hyper_parameters.sigma_w,hyper_parameters.sigma_e,pars,hyper_parameters.center);
end
funv_real_avg=sum(funv_real)/R;
funv_real_noNoise=obj_value_real(hyper_parameters.P,hyper_parameters.Xor,hyper_parameters.B,hyper_parameters.u_real,hyper_parameters.sigma_w,hyper_parameters.sigma_e,pars,hyper_parameters.center); % blurred image before introducing the noise, which is unknown
col={'-o', '-s', '-^', '-v', '-d', '-p'}; % marers for plotting
ticks_to_plot=10; % number of xticks for plotting
legend_plot_f={'SPA','NAM-AG, $s=1$, $r=10$','NAM-AG, $s=10$, $r=10$','NAM-AG, $s=50$, $r=10$','NAM-AG, $s=100$, $r=10$','NAM-AG, $s=200$, $r=10$','$F\left(\mathbf{z}^{\mathrm{real}},\mathbf{u}^{\mathrm{real}}\right)$'}; % methods
legend_plot={'SPA','NAM-AG, $s=1$','NAM-AG, $s=10$','NAM-AG, $s=50$','NAM-AG, $s=100$','NAM-AG, $s=200$'}; % methods

%% Function Values 

funv=zeros(N+1,1+length(AG_iter_s));

for t=1:R
    LoadDir_t=([LoadDir,'\',num2str(t)]);
    
    SPA=load([LoadDir_t,'\SPA.mat']);
    funv(:,1)=funv(:,1)+SPA.SPA.funv/R;
    
    for l=1:length(AG_iter_s)
        s=AG_iter_s(l);
        r=AG_iter_r;
        NAMAG=load([LoadDir_t,'\NAMAG_s',num2str(s),'_r',num2str(r),'.mat']);
        funv(:,l+1)=funv(:,l+1)+NAMAG.NAMAG.funv(1:N+1)/R;
    end
end

%% Plot Function Values (real)

iter_to_plot=25000;
fh=figure(1);

xplot=1:iter_to_plot+1;
for i=1:1+length(AG_iter_s)
    h=semilogy(funv(:,i),col{i}','MarkerIndices',xplot(1:iter_to_plot/ticks_to_plot:end),'LineWidth',2.5);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    set(gca,'FontSize',20)
    hold on
end
%semilogy(funv_real_noNoise*ones(iter_to_plot+1,1),'black','LineWidth',1.5)
hold off

grid on
xlim([0 iter_to_plot+1])
set(gca,'XTickLabel',strsplit(num2str(xticks)))
xlabel('$\mathbf{Total\ Iterations}$','FontSize',24,'Interpreter','latex')
ylabel('$\mathbf{Function\ Value}\ \left(F_j\right)$','FontSize',24,'Interpreter','latex')
%title('$\mathbf{Average\ Function\ Values\ Over\ 100\ Monte\ Carlo\ Trials}$','FontSize',14,'Interpreter','latex')
legend(legend_plot,'FontSize',24,'Interpreter','latex')
hold off
fh.WindowState = 'maximized';
ax = gca;
exportgraphics(ax,'funv_real.eps')

%% Plot Function Values (convegence)

iter_to_plot=25000;
fh=figure(2);

xplot=1:iter_to_plot+1;
for i=1:1+length(AG_iter_s)
    h=semilogy(funv(:,i)-funv(end,i),col{i},'MarkerIndices',xplot(1:iter_to_plot/ticks_to_plot:end),'LineWidth',2.5);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    set(gca,'FontSize',20)
    hold on
end
hold off

grid on
xlim([0 iter_to_plot+1])
set(gca,'XTickLabel',strsplit(num2str(xticks)))
xlabel('$\mathbf{Total\ Iterations}$','FontSize',24,'Interpreter','latex')
ylabel('$\mathrm{Dev}_j^{\mathrm{con}}$','FontSize',24,'Interpreter','latex')
%title('$\mathbf{Average\ Distance\ of\ Function\ Values\ From\ }F\left(\mathbf{z}^{\mathrm{con}},\mathbf{u}^N\right)\mathbf{\ Over\ 100\ Monte\ Carlo\ Trials}$','FontSize',14,'Interpreter','latex')
legend(legend_plot,'FontSize',24,'Interpreter','latex')
hold off
fh.WindowState = 'maximized';
ax = gca;
exportgraphics(ax,'dev.eps')

%% Relative Errors (real)

err_real=zeros(N+1,1+length(AG_iter_s));

for t=1:R
    LoadDir_t=([LoadDir,'\',num2str(t)]);
    
    SPA=load([LoadDir_t,'\SPA.mat']);
    err_real(1:N+1,1)=err_real(1:N+1,1)+SPA.SPA.err(1:N+1)/R;
    
    for l=1:length(AG_iter_s)
        s=AG_iter_s(l);
        r=AG_iter_r;
        NAMAG=load([LoadDir_t,'\NAMAG_s',num2str(s),'_r',num2str(r),'.mat']);
        err_real(1:N+1,l+1)=err_real(1:N+1,l+1)+NAMAG.NAMAG.err(1:N+1)/R;
    end
end

%% Plot Relative Errors (real)

iter_to_plot=25000;
figure(3)

xplot=1:N+1;
for i=1:1+length(AG_iter_s)
    h=plot(err_real(:,i),col{i}','MarkerIndices',xplot(1:iter_to_plot/ticks_to_plot:end),'LineWidth',2.5);
    set(h, 'MarkerFaceColor', get(h,'Color'));
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

err_con=zeros(N+1,1+length(AG_iter_s));

for t=1:R
    LoadDir_t=([LoadDir,'\',num2str(t)]);
    
    SPA=load([LoadDir_t,'\SPA_result.mat']);
    err_con(:,1)=err_con(:,1)+SPA.SPA_result.err/R;
    
    for l=1:length(AG_iter_s)
        s=AG_iter_s(l);
        r=AG_iter_r;
        NAMAG=load([LoadDir_t,'\NAMAG_result_s',num2str(s),'_r',num2str(r),'.mat']);
        err_con(:,l+1)=err_con(:,l+1)+NAMAG.NAMAG_result.err(1:N+1)/R;
    end
end

%% Plot Relative Errors (convergence)

iter_to_plot=25000;
fh=figure(4);

xplot=1:N+1;
for i=1:1+length(AG_iter_s)
    h=plot(err_con(:,i),col{i}','MarkerIndices',xplot(1:iter_to_plot/ticks_to_plot:end),'LineWidth',2.5);
    set(h,'MarkerFaceColor',get(h,'Color'));
    set(gca,'FontSize',20)
    hold on
end
hold off

grid on
xlim([0 iter_to_plot+1])
set(gca,'XTickLabel',strsplit(num2str(xticks)))
xlabel('$\mathbf{Total\ Iterations}$','FontSize',24,'Interpreter','latex')
ylabel('$\mathrm{ConGap}_j$','FontSize',24,'Interpreter','latex')
%title('$\mathbf{Average\ Relative\ Error\ Over\ 100\ Monte\ Carlo\ Trials\ with\ Respect\ to\ z^{con}}$','FontSize',24,'Interpreter','latex')
legend(legend_plot,'FontSize',24,'Interpreter','latex')
hold off
fh.WindowState = 'maximized';
ax = gca;
exportgraphics(ax,'err.eps')

%% Plot and Save Real and Blurred Image of Trial #1

X=double(imread('cameraman.tif'));
plot_image(X/255,5,'$\mathbf{Real\ Image}$')
ax = gca;
exportgraphics(ax,'real_image.eps')

load([LoadDir,'\1\Bobs.mat']);
plot_image(Bobs,6,'$\mathbf{Blurred\ and\ Noisy\ Image}$')
ax = gca;
exportgraphics(ax,'blur_image.eps')

%% Plot Output Images of Trial #1

images_to_plot=1:5;
load([LoadDir,'\1\images.mat']);

fh=figure(7);
mon_SPA=[]; mon_s1=[]; mon_s10=[]; mon_s50=[]; mon_s100=[]; mon_s200=[];
for j=images_to_plot
    mon_SPA=[mon_SPA x_SPA(:,:,j)];
    mon_s1=[mon_s1 x_NAMAG_s1(:,:,j)];
    mon_s10=[mon_s10 x_NAMAG_s10(:,:,j)];
    mon_s50=[mon_s50 x_NAMAG_s50(:,:,j)];
    mon_s100=[mon_s100 x_NAMAG_s100(:,:,j)];
    mon_s200=[mon_s200 x_NAMAG_s200(:,:,j)];
end
mon=[mon_SPA;mon_s1;mon_s10;mon_s50;mon_s100;mon_s200];
montage(mon)
pbaspect([6 5 1])
fh.WindowState = 'maximized';
ax = gca;
exportgraphics(ax,'output_images.eps')
