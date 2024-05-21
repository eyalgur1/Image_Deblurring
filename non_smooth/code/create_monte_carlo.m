%% Generate Random Matrices

R=50; % number of realizations
sigma_w=1e-4; % noising paramter
sigma_e=1e-4; % noising paramter
SaveDirI='C:\Users\eyal.gur\OneDrive - Technion\06 - NSNAM Paper\input\R\';
%Xor=double(imread('camerman.png')); % The camerman picture
Xor=double(imread('gantrycrane.png'))/255; % The gantry crane picture normalized to the interval [0,1]
image(Xor)
Xor=Xor(:,1:256,:); % take only 256x256x3 sub-picture
plot_image(Xor,5,'$\mathbf{Real\ Image}$') % plot the original picture
[P,center]=psfGauss([5,5],2); % Chosing the PSF
[p,S,d]=detect_structure(P); % Detect the structure of the PSF matrix P
B=imfilter(Xor,P,'circular'); % Bluring the given picture X ("real" blurred image)
plot_image(B,6,'$\mathbf{Real\ Image}$') % plot some blurred picture

%%
hyper_parameters=struct;
hyper_parameters.sigma_w=sigma_w;
hyper_parameters.sigma_e=sigma_e;
hyper_parameters.Xor=Xor;
hyper_parameters.P=P;
hyper_parameters.B=B;
hyper_parameters.center=center;
hyper_parameters.u_real=d;

save([SaveDirI,'/hyper_parameters.mat'],'hyper_parameters')

for j=1:R
    SaveDir=[SaveDirI,num2str(j)];
    mkdir(SaveDir);
    
    Pobs=P; % noising the PSF matrix P (differently for any component)
    for i=1:p
        Pobs=Pobs+sigma_e*rand(1)*S(:,:,i);
    end   
    Bobs=B+sigma_w*randn(size(B)); % noising the blurred picture
    
    save([SaveDir,'/Pobs.mat'],'Pobs')
    save([SaveDir,'/Bobs.mat'],'Bobs')
end