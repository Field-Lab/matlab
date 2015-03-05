path2high='S:\data\alexandra\MEA_data\20121023\fit_res_HC\';
path2low='S:\data\alexandra\MEA_data\20121023\fit_res_LC\';
HC=dir([path2high,'*.mat']);
LC=dir([path2low,'*.mat']);
k=min(length(HC),length(LC));
highContr=zeros(k,144);
lowContr=zeros(k,144);
for i=1:k
    name=LC(i).name;
    load([path2high,name]);
    highContr(i,:)=common_res_fit(1,1:144);
    
    load([path2low,name]);
    lowContr(i,:)=common_res_fit(1,1:144);
end

% plot(highContr','k.');
figure
subplot(3,1,1)
hold on
plot(mean(highContr),'r.','MarkerSize',15)
% plot(lowContr','g.')
plot(mean(lowContr),'b.','MarkerSize',15)
legend('high','low')
set(gca,'XTick',6:12:12*12,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*12
    line([i,i],[0, 150],'color','m')    
end
axis([0, 144, 0,150])
title('Latency')



for i=1:length(LC)
    load([path2high,HC(i).name]);
    highContr(i,:)=abs(common_res_fit(1,1:144));
    
    load([path2low,LC(i).name]);
    lowContr(i,:)=abs(common_res_fit(1,1:144));
end
subplot(3,1,2)
hold on
plot(mean(highContr),'r.','MarkerSize',15)
% plot(lowContr','g.')
plot(mean(lowContr),'b.','MarkerSize',15)
legend('high','low')
set(gca,'XTick',6:12:12*12,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*12
    line([i,i],[0, 350],'color','m')    
end
axis([0, 144, 0,350])
title('Amplitude')

for i=1:length(LC)
    load([path2high,HC(i).name]);
    highContr(i,:)=common_res_fit(3,1:144);
    
    load([path2low,LC(i).name]);
    lowContr(i,:)=common_res_fit(3,1:144);
end
subplot(3,1,3)
hold on
plot(mean(highContr),'r.','MarkerSize',15)
% plot(lowContr','g.')
plot(mean(lowContr),'b.','MarkerSize',15)
legend('high','low')
set(gca,'XTick',6:12:12*12,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*12
    line([i,i],[0, 150],'color','m')    
end
axis([0, 144, 0,50])
title('Width')



% plot(highContr','k.');
figure
hold on

path2high='S:\data\alexandra\MEA_data\20121002\fit_res_HC\';
path2low='S:\data\alexandra\MEA_data\20121002\fit_res_LC\';
HC=dir([path2high,'*.mat']);
LC=dir([path2low,'*.mat']);
k=min(length(HC),length(LC));
highContr=zeros(k,144);
lowContr=zeros(k,144);

clear p l t
for i=1:k
    name=LC(i).name;
    load([path2high,name]);
    highContr(i,:)=common_res_fit(1,1:144);
    
    load([path2low,name]);
    lowContr(i,:)=common_res_fit(1,1:144);
end
for i=1:144
    p(1:k,i)=highContr(:,i)~=0&lowContr(:,i)~=0;
end
for i=1:k
    name=LC(i).name;
    load([path2high,name]);
    highContr(i,:)=common_res_fit(2,1:144);
    
    load([path2low,name]);
    lowContr(i,:)=common_res_fit(2,1:144);
end
for i=1:144
    l(i)=mean(highContr(p(:,i),i)-lowContr(p(:,i),i));
    t(i)=std(highContr(p(:,i),i)-lowContr(p(:,i),i))/sqrt(sum(p(:,i)));
end

errorbar(l,t,'b.','MarkerSize',15)
set(gca,'XTick',6:12:12*12,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*12
    line([i,i],[-20, 20],'color','m')    
end
line([0,144],[0, 0],'color','k')    
axis([0, 144, -20,20])

legend('rod','wt')
title('high - low, blue rod, red wt')