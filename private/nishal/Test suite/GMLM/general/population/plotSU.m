function plotSU(K,mask)
K2=gather(K);
mm = reshape(1:3200,[80,40]);
iidx = mm(mask);
[r,c] = find(reshape(mask,[80,40])');


figure;
for isu=1:size(K2,2)
    subplot(7,4,isu);
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