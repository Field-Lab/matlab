clear

totTime=120;
timing=18:40:totTime;
path2fit=['S:\data\alexandra\MEA_data\20120902_1\fit_res_CB\'];
fit=dir([path2fit,'*.mat']);
common=zeros(3,3,36);
for unit=1:length(fit)
    load([path2fit,fit(unit).name]);
    common(:,:,unit)=common_res_fit(:,1:3);    
end
path2fit=['S:\data\alexandra\MEA_data\20120902_2\fit_res_CB\'];
fit=dir([path2fit,'*.mat']);
for unit=1:length(fit)
    load([path2fit,fit(unit).name]);
    common(:,:,unit+21)=common_res_fit;    
end

a=reshape(common(2,:,:),3,36);
b=mean(a,2);
c=std(a,0,2)/sqrt(36);
errorbar(timing,b,c,'.')
hold on

totTime=120;
timing=sort([2:40:totTime 7:40:totTime 23:40:totTime 29:40:totTime]);
path2fit=['S:\data\alexandra\MEA_data\20120902_1\fit_res_HC\'];
fit=dir([path2fit,'*.mat']);
commonFF=zeros(3,12,28);
for unit=1:length(fit)
    load([path2fit,fit(unit).name]);
    commonFF(:,:,unit)=common_res_fit(:,1:12);    
end
path2fit=['S:\data\alexandra\MEA_data\20120902_2\fit_res_HC\'];
fit=dir([path2fit,'*.mat']);
for unit=1:length(fit)
    load([path2fit,fit(unit).name]);
    commonFF(:,:,unit+18)=common_res_fit;    
end

a=reshape(commonFF(2,:,:),12,28);
b=mean(a,2);
c=std(a,0,2)/sqrt(28);
errorbar(timing,b,c,'.r')

for i=40:40:totTime
    line([i,i],[40,140],'color','k')
end
axis([0 totTime 60 120])

set(gca,'xtick',20:40:totTime,'xticklabel',{'4','3','2'})
xlabel('ND')
ylabel('latency, ms')
title([{'WT latency, mean +-s.e.m, 2 retinas'},{'red - full field (28 units), blue - checkerboard (36 units)'}])