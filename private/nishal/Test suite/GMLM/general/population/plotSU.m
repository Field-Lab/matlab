function h=plotSU(K,mask)
K2=gather(K);
mm = reshape(1:3200,[80,40]);
iidx = mm(mask);
[r,c] = find(reshape(mask,[80,40])');

nSU=size(K2,2);
h=figure;
for isu=1:nSU
    subplot(ceil(sqrt(nSU)),ceil(sqrt(nSU)),isu);
    u=zeros(3200,1);
    u(iidx) = K2(:,isu);
    xx = reshape(u,[80,40])';
    
    imagesc(xx(min(r):max(r),min(c):max(c)));
    axis image
    set(gca,'xTick',[]);
    set(gca,'yTick',[]);
    colormap gray
end

end