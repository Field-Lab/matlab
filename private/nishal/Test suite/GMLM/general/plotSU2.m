function plotSU2(K,mask,dim1,dim2)
K2=gather(K);
mm = reshape(1:dim1*dim2,[dim1,dim2]);
iidx = mm(mask);
[r,c] = find(reshape(mask,[dim1,dim2])');


figure;
for isu=1:size(K2,2)
    subplot(7,4,isu);
    u=zeros(dim1*dim2,1);
    u(iidx) = K2(:,isu);
    xx = reshape(u,[dim1,dim2])';
    
    imagesc(xx(min(r):max(r),min(c):max(c)));
    axis image
    set(gca,'xTick',[]);
    set(gca,'yTick',[]);
    colormap gray
end

end