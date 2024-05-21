%% Generate Random Matrices

R=100; % number of realizations
sigma_w=1e-4; % noising paramter
sigma_e=1e-3; % noising paramter
SaveDirO='C:/Users/eyal.gur/OneDrive - Technion/04 - NAM Paper/14 - Numerical Analysis - 27.6.21/R/';
Xor=double(imread('cameraman.tif')); % The camerman picture
Xor=Xor/255; % Normalization of the picture to the interval [0,1]
[P,center]=psfGauss([5,5],2); % Chosing the PSF
[p,S,d]=detect_structure(P); % Detect the structure of the PSF matrix P
B=imfilter(Xor,P,'circular'); % Bluring the given picture X

%%
hyper_parameters=struct;
hyper_parameters.sigma_w=sigma_w;
hyper_parameters.sigma_e=sigma_e;
hyper_parameters.Xor=Xor;
hyper_parameters.P=P;
hyper_parameters.B=B;
hyper_parameters.center=center;
hyper_parameters.u_real=d;

save([SaveDirO,'/hyper_parameters.mat'],'hyper_parameters')

for j=1:R
    SaveDir=[SaveDirO,num2str(j)];
    mkdir(SaveDir);
    
    Pobs=P; % noising the PSF matrix P (differently for any component)
    for i=1:p
        Pobs=Pobs+sigma_e*rand(1)*S(:,:,i);
    end   
    Bobs=B+sigma_w*randn(size(B)); % noising the blurred picture
    
    save([SaveDir,'/Pobs.mat'],'Pobs')
    save([SaveDir,'/Bobs.mat'],'Bobs')
end