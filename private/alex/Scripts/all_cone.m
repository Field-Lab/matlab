clear
path2coef='S:\user\alexandra\MEA_data\20120920\fit_res_HC\';
nd='654321234'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    cone(i,:,1:3)=[zeros(24,3); common_res_fit'; zeros(21,3)];
end

clear res
path2coef='S:\user\alexandra\MEA_data\20120920_1\fit_res_HC\';
nd='543212'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=[zeros(36,3); common_res_fit'; zeros(45,3)];
end
cone=[cone; res];

clear res
path2coef='S:\user\alexandra\MEA_data\20120921\fit_res_HC\';
nd='543212345'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=[zeros(36,3); common_res_fit(:,1:89)'; zeros(21,3)];
end
cone=[cone; res];

clear res
path2coef='S:\user\alexandra\MEA_data\20120925\fit_res_HC\';
nd='654323232123'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=[zeros(24,3); common_res_fit(:,1:60)';zeros(62,3)];
end
cone=[cone; res];




clear rod
path2coef='S:\user\alexandra\MEA_data\20120928\fit_res_HC\';
nd='7654321234'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    rod(i,:,1:3)=[zeros(12,3); common_res_fit'; zeros(22,3)];
end

clear res
path2coef='S:\user\alexandra\MEA_data\20121002\fit_res_HC\';
nd='876543212345'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end
rod=[rod; res];


clear wt
path2coef='S:\user\alexandra\MEA_data\20121023\fit_res_HC\';
nd='87654321234567'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    wt(i,:,1:3)=common_res_fit(:,1:146)';
end

clear res common_res_fit coef path2coef






figure
set(gcf,'position',[1 31 1680 946])

subplot(3,1,1)
plot(wt(:,:,2)','b.')
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[0, 200],'color','m')    
end
a=wt(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'k.','MarkerSize',20)
axis([0,147,0,200])
title('Latency')

subplot(3,1,2)
plot(wt(:,:,1)','b.')
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[-700, 700],'color','m')    
end
a=wt(:,:,1)';
for i=1:size(a,1)
    tmp= a(i,:)>0;
    b(i)=mean(a(i,tmp));
    tmp=find(a(i,:)<0&a(i,:)>-500);
    c(i)=mean(a(i,tmp));
end
plot(b,'k.','MarkerSize',20)
plot(c,'k.','MarkerSize',20)
line([0,147],[0,0],'color','k')
axis([0,147,-700,700])
title('Amplitude')

subplot(3,1,3)
plot(wt(:,:,3)','b.')
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[0, 80],'color','m')    
end
a=wt(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'k.','MarkerSize',20)
line([0,147],[0,0],'color','k')
axis([0,147,0,80])
title('Width')






figure
set(gcf,'position',[1 31 1680 946])
a=sum(wt(:,:,1),2);
onCell=a>0;
offCell=a<0;
a=sum(cone(:,:,1),2);
onCellcone=a>0;
offCellcone=a<0;
a=sum(rod(:,:,1),2);
onCellrod=a>0;
offCellrod=a<0;

subplot(3,1,1)
hold on
plot(wt(onCell,:,2)','b.')
plot(wt(offCell,:,2)','r.')
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[0, 200],'color','m')    
end
a=wt(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'k.','MarkerSize',20)

a=cone(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'g.','MarkerSize',20)

a=rod(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'m.','MarkerSize',20)

axis([0,147,0,200])
title('Latency')

subplot(3,1,2)
hold on
plot(wt(onCell,:,1)','b.')
plot(wt(offCell,:,1)','r.')
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[-700, 700],'color','m')    
end
a=wt(:,:,1)';
for i=1:size(a,1)
    tmp= a(i,:)>0;
    b(i)=mean(a(i,tmp));
    tmp=find(a(i,:)<0&a(i,:)>-500);
    c(i)=mean(a(i,tmp));
end
plot(b,'k.','MarkerSize',20)
plot(c,'k.','MarkerSize',20)

a=cone(:,:,1)';
for i=1:size(a,1)
    tmp= a(i,:)>0;
    b(i)=mean(a(i,tmp));
    tmp=find(a(i,:)<0&a(i,:)>-500);
    c(i)=mean(a(i,tmp));
end
plot(b,'g.','MarkerSize',20)
plot(c,'g.','MarkerSize',20)

a=rod(:,:,1)';
for i=1:size(a,1)
    tmp= a(i,:)>0;
    b(i)=mean(a(i,tmp));
    tmp=find(a(i,:)<0&a(i,:)>-500);
    c(i)=mean(a(i,tmp));
end
plot(b,'m.','MarkerSize',20)
plot(c,'m.','MarkerSize',20)

line([0,147],[0,0],'color','k')
axis([0,147,-700,700])
title('Amplitude')

subplot(3,1,3)
hold on
plot(wt(onCell,:,3)','b.')
plot(wt(offCell,:,3)','r.')
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[0, 80],'color','m')    
end
a=wt(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
plot(b,'k.','MarkerSize',20)

a=cone(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
plot(b,'g.','MarkerSize',20)

a=rod(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
plot(b,'m.','MarkerSize',20)

line([0,147],[0,0],'color','k')
axis([0,147,0,80])
title('Width')














figure
set(gcf,'position',[1 31 1680 946])

subplot(3,1,1)
plot(wt(:,:,2)','b.')
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[0, 200],'color','m')    
end
a=wt(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'k.','MarkerSize',20)
axis([0,147,0,200])
title('Latency')

subplot(3,1,2)
plot(wt(:,:,1)','b.')
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[-700, 700],'color','m')    
end
a=wt(:,:,1)';
for i=1:size(a,1)
    tmp= a(i,:)>0;
    b(i)=mean(a(i,tmp));
    tmp=find(a(i,:)<0&a(i,:)>-500);
    c(i)=mean(a(i,tmp));
end
plot(b,'k.','MarkerSize',20)
plot(c,'k.','MarkerSize',20)
line([0,147],[0,0],'color','k')
axis([0,147,-700,700])
title('Amplitude')

subplot(3,1,3)
plot(wt(:,:,3)','b.')
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[0, 80],'color','m')    
end
a=wt(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'k.','MarkerSize',20)
line([0,147],[0,0],'color','k')
axis([0,147,0,80])
title('Width')






figure
set(gcf,'position',[1 31 1680 946])
a=sum(wt(:,:,1),2);
onCell=a>0;
offCell=a<0;
a=sum(cone(:,:,1),2);
onCellcone=a>0;
offCellcone=a<0;
a=sum(rod(:,:,1),2);
onCellrod=a>0;
offCellrod=a<0;

subplot(3,1,1)
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[0, 200],'color','m')    
end
a=wt(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'k.','MarkerSize',20)

a=cone(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'g.','MarkerSize',20)

a=rod(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'m.','MarkerSize',20)

axis([0,147,0,200])
title('Latency')

subplot(3,1,2)
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[-700, 700],'color','m')    
end
a=wt(:,:,1)';
for i=1:size(a,1)
    tmp= a(i,:)>0;
    b(i)=mean(a(i,tmp));
    tmp=find(a(i,:)<0&a(i,:)>-500);
    c(i)=mean(a(i,tmp));
end
plot(b,'k.','MarkerSize',20)
plot(c,'k.','MarkerSize',20)

a=cone(:,:,1)';
for i=1:size(a,1)
    tmp= a(i,:)>0;
    b(i)=mean(a(i,tmp));
    tmp=find(a(i,:)<0&a(i,:)>-500);
    c(i)=mean(a(i,tmp));
end
plot(b,'g.','MarkerSize',20)
plot(c,'g.','MarkerSize',20)

a=rod(:,:,1)';
for i=1:size(a,1)
    tmp= a(i,:)>0;
    b(i)=mean(a(i,tmp));
    tmp=find(a(i,:)<0&a(i,:)>-500);
    c(i)=mean(a(i,tmp));
end
plot(b,'m.','MarkerSize',20)
plot(c,'m.','MarkerSize',20)

line([0,147],[0,0],'color','k')
axis([0,147,-700,700])
title('Amplitude')

subplot(3,1,3)
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[0, 80],'color','m')    
end
a=wt(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
plot(b,'k.','MarkerSize',20)

a=cone(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
plot(b,'g.','MarkerSize',20)

a=rod(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
plot(b,'m.','MarkerSize',20)

line([0,147],[0,0],'color','k')
axis([0,147,0,80])
title('Width')





figure
set(gcf,'position',[1 31 1680 946])
tmp=wt(:,12:12:end,2)';
Lat=diff(tmp);
tmp=abs(wt(:,12:12:end,1)');
Amp=diff(tmp);

tmp=cone(:,12:12:end,2)';
Lat1=diff(tmp);
tmp=abs(cone(:,12:12:end,1)');
Amp1=diff(tmp);

tmp=rod(:,12:12:end,2)';
Lat2=diff(tmp);
tmp=abs(rod(:,12:12:end,1)');
Amp2=diff(tmp);

nd='87654321234567'
for i=1:9
    subplot(3,3,i)
    hold on    
    plot(Lat(i,:),Amp(i,:),'.b')
    plot(Lat1(i,:),Amp1(i,:),'.g')
    plot(Lat2(i,:),Amp2(i,:),'.r')
    if i==1
        legend('wt','cone','rod')
    end
    title(['ND',nd(i+1),' - ND',nd(i),' difference of last trials in ND'])
%     axis([-1.5 1.5 -2 2])
    line([-150 150],[0 0],'color','k')
    line([0 0],[-300 300],'color','k')
    xlabel('Latency')
    ylabel('Amplitude')
    axis([-150 150 -300 300])
end

