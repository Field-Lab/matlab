function plot_SU_with_idx(sta,mask)

[r,c] =find(reshape(mask,[40,40])');
% figure;
imagesc(((sta)));
colormap gray
colorbar
axis image
hold on;
iidx = reshape(1:1600,[40,40]);
xx = repmat([1:40]',[1,40]); xx=xx(mask);
yy=repmat(1:40,[40,1]); yy =yy(mask);
pixidx =iidx(mask);
for ipix=1:length(pixidx)
    hold on;
text(xx(ipix)-0.5,yy(ipix),sprintf('%d',(ipix)),'Color','r','FontSize',8);
end