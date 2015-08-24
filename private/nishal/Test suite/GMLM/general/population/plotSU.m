function plotSU(K,mask)
K2=gather(K);
mm = reshape(1:3200,[80,40]);
iidx = mm(mask);


figure;
for isu=1:size(K2,2)
    subplot(5,2,isu);
    u=zeros(3200,1);
    u(iidx) = K2(:,isu);
    imagesc(reshape(u,[80,40])');
    axis image
    set(gca,'xTick',[]);
    set(gca,'yTick',[]);
end

end