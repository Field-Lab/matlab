clear
path2fit=['S:\data\alexandra\MEA_data\20120830_1\fit_res_CBF\'];
fit=dir([path2fit,'*.mat']);

nds='765432'
common=zeros(3,16,10);
for unit=1:length(fit)
    load([path2fit,fit(unit).name]);
    common(:,:,unit)=common_res_fit;    
end
a=reshape(common(2,:,:),16,10);
b=mean(a,2);
c=std(a,0,2)/sqrt(10);
errorbar(m,b,c)

m=sort([4:28:168 12:28:150 20:28:150]);
figure
plot(a,'.')

for i=1:16
    tmp=a(i,:)>0;
    allMean(i)=mean(a(i,tmp));
    allSem(i)=std(a(i,tmp)/sqrt(sum(tmp)));
end
figure
errorbar(m,allMean,allSem)
set(gca,'xtick', 14:28:168,'xticklabel',{'ND7','ND6','ND5','ND4','ND3','ND2'})
title({'Melanopsin checkerboard temporal development','1 retina, 10 units'})