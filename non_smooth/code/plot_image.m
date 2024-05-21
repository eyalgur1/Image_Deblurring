% Plots a 256X256 image with the correct scale, color and labels.
% X - the required image to plot (256x256 matrix)
% index - the index of the image
% cmap - colormap of the photo (obtained by the command colormap given after any plot) 
% tit - a script for the title of the image

function plot_image(X,index,tit)
figure(index)
imshow(X)
%colormap(cmap)
%yticks(1:50:256)
%yticklabels({'250','200','150','100','50','0'})
axis square
xlim([1 256])
ylim([1 256])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title(tit,'interpreter','latex','FontSize',14)
hold off
end