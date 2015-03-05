
figure
clear
path2coef='S:\user\alexandra\MEA_data\20120329\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
res=zeros(length(coef),144,3);
fillIn=sort([2:12:144 4:12:144 10:12:140 12:12:140]);
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,fillIn,1:3)=common_res_fit';
end

subplot(3,1,1)
plot(res(:,:,2)','b.')
set(gca,'XTick',6:12:12*12,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*12
    line([i,i],[0, 200],'color','m')    
end
a=res(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'k.','MarkerSize',20)
axis([0,144,0,200])
title('Latency')

subplot(3,1,2)
plot(res(:,:,1)','b.')
set(gca,'XTick',6:12:12*12,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*12
    line([i,i],[-500, 500],'color','m')    
end
line([0,144],[0,0],'color','k')
axis([0,144,-500,500])
title('Amplitude')

subplot(3,1,3)
plot(res(:,:,3)','b.')
set(gca,'XTick',6:12:12*12,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*12
    line([i,i],[0, 80],'color','m')    
end
line([0,144],[0,0],'color','k')
axis([0,144,0,80])
title('Width')