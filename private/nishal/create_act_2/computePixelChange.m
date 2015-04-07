function [proj,maxChange,averageChange]=computePixelChange(mov_new,mov_orig,A,Ainv,CellMasks_sp)

movLen=size(mov_orig,3);

proj=zeros(movLen,1);
maxChange=zeros(movLen,1);
averageChange=zeros(movLen,1);
for iframe=1:movLen
    mov_fr = mov_orig(:,:,iframe);
  
    
    mov_new_fr = mov_new(:,:,iframe);
  
    proj(iframe)=norm(Ainv*A*mov_fr(:));
    mov_fr=mov_fr(logical(CellMasks_sp{1}));
    mov_new_fr=mov_new_fr(logical(CellMasks_sp{1}));
    
    maxChange(iframe)=max(max(abs(mov_fr-mov_new_fr)));
    averageChange(iframe)=mean(mean(abs(mov_new_fr-mov_fr)));
end

idx=1:movLen;
valididx=idx>0.1*movLen & idx<0.9*movLen;
proj=proj(valididx);
maxChange=maxChange(valididx);
averageChange=averageChange(valididx);

figure('Color','w');
subplot(3,1,1);
plot(proj,maxChange,'.','MarkerSize',5);
xlabel('Original Projection');
ylabel('Max Change');

subplot(3,1,2);
plot(proj,averageChange,'.','MarkerSize',5);
xlabel('Original Projection');
ylabel('Average Change');

subplot(3,1,3);
hist(proj,100);
title(sprintf('Proj. histogram , percent of fr. with avg change<0.5= %.2f and max change <0.5 = %.2f ',100*sum(averageChange<0.5)/numel(averageChange),100*sum(maxChange<0.5)/numel(maxChange)));

end