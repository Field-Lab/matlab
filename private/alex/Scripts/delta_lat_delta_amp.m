clear
path2coef='S:\user\alexandra\MEA_data\20121023\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end
figure
set(gcf,'position',[1 31 1680 946])
tmp=res(:,12:12:end,2)';
tmp=tmp./(repmat(tmp(3,:),14,1));
Lat=diff(tmp);
tmp=abs(res(:,12:12:end,1)');
tmp=tmp./(repmat(tmp(3,:),14,1));
Amp=diff(tmp);
nd='87654321234567'
for i=1:9
    subplot(3,3,i)
    hold on
    plot(Lat(i,:),Amp(i,:),'.b')
    title(['ND',nd(i+1),' - ND',nd(i),' difference of last trials in ND'])
%     axis([-1.5 1.5 -2 2])
    line([-1.5 1.5],[0 0],'color','k')
    line([0 0],[-2 2],'color','k')
    xlabel('Latency')
    ylabel('Amplitude')
end


path2coef='S:\user\alexandra\MEA_data\20121002\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
clear res
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end
nd='876543212345'

tmp=res(:,12:12:end,2)';
tmp=tmp./(repmat(tmp(3,:),12,1));
Lat=diff(tmp);
tmp=abs(res(:,12:12:end,1)');
tmp=tmp./(repmat(tmp(3,:),12,1));
Amp=diff(tmp);
nd='87654321234567'
for i=1:9
    subplot(3,3,i)
    plot(Lat(i,:),Amp(i,:),'.r')
end




path2coef='S:\user\alexandra\MEA_data\20120920_1\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
clear res
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end
nd='54321'

tmp=res(:,12:12:end,2)';
tmp=tmp./(repmat(tmp(3,:),5,1));
Lat=diff(tmp);
tmp=abs(res(:,12:12:end,1)');
tmp=tmp./(repmat(tmp(3,:),5,1));
Amp=diff(tmp);
for i=4:7
    subplot(3,3,i)
    plot(Lat(i-3,:),Amp(i-3,:),'.g')
end








clear
path2coef='S:\user\alexandra\MEA_data\20121023\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end

subplot(3,1,1)
plot(res(:,:,2)','b.','MarkerSize',1)
set(gca,'XTick',6:12:12*14,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*14
    line([i,i],[0, 200],'color','m')    
end
a=res(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    d(i)=mean(a(i,tmp));
end
hold on
plot(d,'k.','MarkerSize',20)
line([0,160],[0,0],'color','k')
axis([0,160,0,200])
title('Latency')

subplot(3,1,2)
plot(res(:,:,1)','b.','MarkerSize',1)
set(gca,'XTick',6:12:12*14,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*14
    line([i,i],[-700, 700],'color','m')    
end
line([0,160],[0,0],'color','k')
axis([0,160,-700,700])
hold on
a=res(:,:,1)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    d(i)=mean(abs(a(i,tmp)));
end
plot(d,'k.','MarkerSize',20)
title('Amplitude')

subplot(3,1,3)
plot(res(:,:,3)','b.','MarkerSize',1)
set(gca,'XTick',6:12:12*14,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*14
    line([i,i],[0, 80],'color','m')    
end
line([0,160],[0,0],'color','k')
axis([0,160,0,80])
hold on
a=res(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    d(i)=mean(a(i,tmp));
end
plot(d,'k.','MarkerSize',20)
title('Width')


figure

clear res
path2coef='S:\user\alexandra\MEA_data\20121002\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end

subplot(3,1,1)
plot(res(:,:,2)','m.','MarkerSize',5)

a=res(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
% plot(b,'r.','MarkerSize',20)

subplot(3,1,2)
plot(res(:,:,1)','m.','MarkerSize',5)
a=res(:,:,1)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(abs(a(i,tmp)));
end
% plot(b,'r.','MarkerSize',20)


subplot(3,1,3)
plot(res(:,:,3)','m.','MarkerSize',5)
a=res(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
% plot(b,'r.','MarkerSize',20)


clear res
path2coef='S:\user\alexandra\MEA_data\20120920_1\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end

subplot(3,1,1)
plot(37:length(res)+36,res(:,:,2)','g.','MarkerSize',1)
a=res(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(a(i,tmp));
end
hold on
plot(37:length(res)+36,c,'c.','MarkerSize',20)

subplot(3,1,2)
plot(37:length(res)+36,res(:,:,1)','g.','MarkerSize',1)
a=res(:,:,1)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(abs(a(i,tmp)));
end
plot(37:length(res)+36,c,'c.','MarkerSize',20)

subplot(3,1,3)
plot(37:length(res)+36,res(:,:,3)','g.','MarkerSize',1)
a=res(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(a(i,tmp));
end
plot(37:length(res)+36,c,'c.','MarkerSize',20)


