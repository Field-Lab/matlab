function plot_SU_with_idx2(sta,mask2)


d1 = size(sta,1);d2=size(sta,2);
[r,c] =find((mask2));

% figure;
imagesc(((sta)));
colormap gray
colorbar
axis image
hold on;
iidx = reshape(1:d1*d2,[d2,d1]);
yy = repmat([1:d1],[d2,1]); yy=yy(mask2);
xx=repmat([1:d2]',[1,d1]); xx =xx(mask2);
pixidx =iidx(mask2);
for ipix=1:length(pixidx)
    hold on;
text(xx(ipix)-0.5,yy(ipix),sprintf('%d',(ipix)),'Color','r','FontSize',8);
end