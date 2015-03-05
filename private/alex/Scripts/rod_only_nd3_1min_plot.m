clear
date='20121004'

path2fit=['S:\data\alexandra\MEA_data\',date,'\fit_res_HC\'];
fits=dir([path2fit,'\*.mat']);
pulled=[];
for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    name=fits(unit).name(1:end-15);
    cellType=sum(sign(common_res_fit(1,:)));
    pulled=[pulled; common_res_fit];  
end
k=pulled(1:3:end,:);
a=sum(k,2);
a=a>0
amplOn=mean(k(a,:));
amplOff=mean(k(~a,:));
lat=mean(pulled(2:3:end,:));
wid=mean(pulled(3:3:end,:));

t=[12 30:9:71 71 72:12:194];
tot_high=[25 9 9 9 9 9 1 12 12 12 12 12 9 9 9 9 9 9 9];
a=cumsum(tot_high);
b(1)=12.5;
for i=1:18
    b(i+1)=a(i)+tot_high(i+1)/2;
end
subplot(3,1,1)
hold on
set(gca,'XTick',b,'xticklabel',{'9','8','7','6','5','4','3','4','3','2','1','2','3','4','5','6','7','8','3'})

plot(lat,'.','MarkerSize',10)
t1=a+0.5;
for i=t1
    line([i,i],[0, 200],'color','k')
end
axis([0,194,0,200])
title('Latency')


subplot(3,1,2)
hold on
set(gca,'XTick',b,'xticklabel',{'9','8','7','6','5','4','3','4','3','2','1','2','3','4','5','6','7','8','3'})

plot(amplOn,'.','MarkerSize',10)
hold on
plot(amplOff,'.','MarkerSize',10)

t1=a+0.5;
for i=t1
    line([i,i],[-250, 250],'color','k')
end
line([0,194],[0,0],'color','k')
axis([0,194,-250,250])
title('Amplitude')


subplot(3,1,3)
hold on
set(gca,'XTick',b,'xticklabel',{'9','8','7','6','5','4','3','4','3','2','1','2','3','4','5','6','7','8','3'})

plot(wid,'.','MarkerSize',10)
t1=a+0.5;
for i=t1
    line([i,i],[0, 50],'color','k')
end
axis([0,194,0,50])
title('Width')

