figure
clear res
path2coef='S:\user\alexandra\MEA_data\20120920\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end

subplot(3,1,1)
hold on
plot(res(:,:,2)','b.','MarkerSize',1)
a=res(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(a(i,tmp));
end
plot(c,'c.','MarkerSize',20)

subplot(3,1,2)
hold on
plot(res(:,:,1)','b.','MarkerSize',1)
a=res(:,:,1)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(abs(a(i,tmp)));
end
plot(c,'c.','MarkerSize',20)

subplot(3,1,3)
hold on
plot(res(:,:,3)','b.','MarkerSize',1)
a=res(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(a(i,tmp));
end
plot(c,'c.','MarkerSize',20)




clear res b c d
path2coef='S:\user\alexandra\MEA_data\20120920_1\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end

subplot(3,1,1)
hold on
plot(13:length(res)+12,res(:,:,2)','r.','MarkerSize',1)
a=res(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(a(i,tmp));
end
plot(13:length(res)+12,c,'m.','MarkerSize',20)

subplot(3,1,2)
hold on
plot(13:length(res)+12,res(:,:,1)','r.','MarkerSize',1)
a=res(:,:,1)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(abs(a(i,tmp)));
end
plot(13:length(res)+12,c,'m.','MarkerSize',20)

subplot(3,1,3)
hold on
plot(13:length(res)+12,res(:,:,3)','r.','MarkerSize',1)
a=res(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(a(i,tmp));
end
plot(13:length(res)+12,c,'m.','MarkerSize',20)




clear res b c d
path2coef='S:\user\alexandra\MEA_data\20120921\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end

subplot(3,1,1)
hold on
plot(13:length(res)+12,res(:,:,2)','g.','MarkerSize',1)
a=res(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(a(i,tmp));
end
plot(13:length(res)+12,c,'g.','MarkerSize',20)

subplot(3,1,2)
hold on
plot(13:length(res)+12,res(:,:,1)','g.','MarkerSize',1)
a=res(:,:,1)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(abs(a(i,tmp)));
end
plot(13:length(res)+12,c,'g.','MarkerSize',20)

subplot(3,1,3)
hold on
plot(13:length(res)+12,res(:,:,3)','g.','MarkerSize',1)
a=res(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(a(i,tmp));
end
plot(13:length(res)+12,c,'g.','MarkerSize',20)



clear res b c d
path2coef='S:\user\alexandra\MEA_data\20120925\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end
res=res(:,1:60,:);

subplot(3,1,1)
hold on
plot(res(:,:,2)','k.','MarkerSize',1)
a=res(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(a(i,tmp));
end
plot(c,'k.','MarkerSize',20)
axis([0 97 0 200])
set(gca,'XTick',6:12:12*8,'xticklabel',{'6','5','4','3','2','1','2'})
for i=12.5:12:12*8
    line([i,i],[0, 200],'color','k')    
end
line([0,97],[0,0],'color','k')



subplot(3,1,2)
hold on
plot(res(:,:,1)','k.','MarkerSize',1)
a=res(:,:,1)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(abs(a(i,tmp)));
end
plot(c,'k.','MarkerSize',20)
axis([0 97 -400 500 ])
set(gca,'XTick',6:12:12*8,'xticklabel',{'6','5','4','3','2','1','2'})
for i=12.5:12:12*8
    line([i,i],[-400, 500],'color','k')    
end
line([0,97],[0,0],'color','k')



subplot(3,1,3)
hold on
plot(res(:,:,3)','k.','MarkerSize',1)
a=res(:,:,3)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    c(i)=mean(a(i,tmp));
end
plot(c,'k.','MarkerSize',20)
axis([ 0 97 0 75])
set(gca,'XTick',6:12:12*8,'xticklabel',{'6','5','4','3','2','1','2'})
for i=12.5:12:12*8
    line([i,i],[0, 75],'color','k')    
end
line([0,97],[0,0],'color','k')
