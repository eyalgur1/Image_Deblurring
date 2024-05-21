%%
s=10;
r=10;
N=3000;
LoadDir='C:\Users\eyal.gur\OneDrive - Technion\04 - NAM Paper\14 - Numerical Analysis - 27.6.21\R';
load([LoadDir,'\hyper_parameters.mat']);
load([LoadDir,'\1\Bobs.mat']);
load([LoadDir,'\1\Pobs.mat']);
reg_g=@(X)norm(X,2)^2;
clear pars;
pars.regfun=reg_g;
pars.regtran=reg_L;
pars.lambda=0.02;
pars.BC='periodic';

x=NAMAG_RSTLS_video(Pobs,Bobs,hyper_parameters.center,hyper_parameters.sigma_w,...
    hyper_parameters.sigma_e,N,pars,s,r);

%%
v = VideoWriter('u','Motion JPEG AVI');
v.FrameRate=150;
v.Quality=100;
open(v);
for k = 1:150
    k
    imagesc(x(:,:,1))
    colormap(gray)
    axis equal
    xlim([0 256])
    ylim([0 256])
    frame = getframe;
    writeVideo(v,frame);
end
for k = 1:300
    k
    imagesc(x(:,:,k))
    colormap(gray)
    axis equal
    xlim([0 256])
    ylim([0 256])
    frame = getframe;
    writeVideo(v,frame);
end
for k = 300:30:3000
    k
    imagesc(x(:,:,k))
    colormap(gray)
    axis equal
    xlim([0 256])
    ylim([0 256])
    frame = getframe;
    writeVideo(v,frame);
end
for k = 1:150
    k
    imagesc(x(:,:,N))
    colormap(gray)
    axis equal
    xlim([0 256])
    ylim([0 256])
    frame = getframe;
    writeVideo(v,frame);
end
close(v);

% Give the first frames more time
% do more than 500 frames
% set better ratio

%%
surf(peaks)
cdata = print('-RGBImage');
imshow(cdata)