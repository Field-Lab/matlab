clear
path2coef='S:\data\alexandra\MEA_data\20120920\fit_res_HC\';
nd='654321234'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    cone(i,:,1:3)=[zeros(24,3); common_res_fit'; zeros(21,3)];
end

clear res
path2coef='S:\data\alexandra\MEA_data\20120920_1\fit_res_HC\';
nd='543212'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=[zeros(36,3); common_res_fit'; zeros(45,3)];
end
cone=[cone; res];

clear res
path2coef='S:\data\alexandra\MEA_data\20120921\fit_res_HC\';
nd='543212345'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=[zeros(36,3); common_res_fit(:,1:89)'; zeros(21,3)];
end
cone=[cone; res];

clear res
path2coef='S:\data\alexandra\MEA_data\20120925\fit_res_HC\';
nd='654323232123'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=[zeros(24,3); common_res_fit(:,1:60)';zeros(62,3)];
end
cone=[cone; res];




clear rod
path2coef='S:\data\alexandra\MEA_data\20120928\fit_res_HC\';
nd='7654321234'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    rod(i,:,1:3)=[zeros(12,3); common_res_fit'; zeros(22,3)];
end

clear res
path2coef='S:\data\alexandra\MEA_data\20121002\fit_res_HC\';
nd='876543212345'
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end
rod=[rod; res];


clear wt
path2coef='S:\data\alexandra\MEA_data\20121023\fit_res_HC\';
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

a=wt(:,:,2)';
b=cone(:,:,2)';
c=rod(:,:,2)';
clear k
for i=1:size(a,1)
    tmp=a(i,:)>0&onCell';
    k(i,1)=mean(a(i,tmp));
    tmp=a(i,:)>0&offCell';
    k(i,2)=mean(a(i,tmp));
    tmp=b(i,:)>0&onCellcone';
    k(i,3)=mean(b(i,tmp));
    tmp=b(i,:)>0&offCellcone';
    k(i,4)=mean(b(i,tmp));
    tmp=c(i,:)>0&onCellrod';
    k(i,5)=mean(c(i,tmp));
    tmp=c(i,:)>0&offCellrod';
    k(i,6)=mean(c(i,tmp));
end
k(isnan(k))=0;

plot(k,'.','MarkerSize',10)
legend({'wt ON 26', 'wt OFF 9','cone ON 52','cone OFF 39', 'rod ON 20','rod OFF 33'})

for i=12.5:12:146
    line([i,i],[0, 200],'color','k')    
end
axis([0,147,50,200])
title('Latency')

subplot(3,1,2)
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[-400, 400],'color','k')    
end


a=wt(:,:,1)';
b=cone(:,:,1)';
c=rod(:,:,1)';
for i=1:size(a,1)
    tmp=a(i,:)>0&onCell';
    m(i,1)=mean(a(i,tmp));
    tmp=a(i,:)<0&offCell'&a(i,:)>-500;
    m(i,2)=mean(a(i,tmp));
    tmp=b(i,:)>0&onCellcone';
    m(i,3)=mean(b(i,tmp));
    tmp=b(i,:)<0&offCellcone';
    m(i,4)=mean(b(i,tmp));
    tmp=c(i,:)>0&onCellrod';
    m(i,5)=mean(c(i,tmp));
    tmp=c(i,:)<0&offCellrod';
    m(i,6)=mean(c(i,tmp));
end
m(isnan(m))=0;

plot(m,'.','MarkerSize',10)

line([0,147],[0,0],'color','k')
axis([0,147,-400,400])
title('Amplitude')

subplot(3,1,3)
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[0, 80],'color','k')    
end

a=wt(:,:,3)';
b=cone(:,:,3)';
c=rod(:,:,3)';
for i=1:size(a,1)
    tmp=a(i,:)>0&onCell';
    n(i,1)=mean(a(i,tmp));
    tmp=a(i,:)>0&offCell';
    n(i,2)=mean(a(i,tmp));
    tmp=b(i,:)>0&onCellcone';
    n(i,3)=mean(b(i,tmp));
    tmp=b(i,:)>0&offCellcone';
    n(i,4)=mean(b(i,tmp));
    tmp=c(i,:)>0&onCellrod';
    n(i,5)=mean(c(i,tmp));
    tmp=c(i,:)>0&offCellrod';
    n(i,6)=mean(c(i,tmp));
end
n(isnan(n))=0;

plot(n,'.','MarkerSize',10)
line([0,147],[0,0],'color','k')
axis([0,147,0,80])
title('Width')





a=wt(onCell,:,2)';
b=cone(onCellcone,:,2)';
c=rod(onCellrod,:,2)';
a1=wt(onCell,:,1)';
b1=cone(onCellcone,:,1)';
c1=rod(onCellrod,:,1)';
a=wt(:,:,2)';
b=cone(:,:,2)';
c=rod(:,:,2)';
a1=abs(wt(:,:,1)');
b1=abs(cone(:,:,1)');
c1=abs(rod(:,:,1)');

cnt=1;
clear all_corr
for i=1:12:144
    tmp=[];
    cnt1=1;
    for j=1:size(a,2)
        tt=corr(a(i:i+11,j),a1(i:i+11,j));
        if ~isnan(tt)
            tmp(cnt1)=tt;
            cnt1=cnt1+1;
        end
    end
    if ~isempty(tmp)
        all_corr(1,cnt)=mean(tmp);
    end
    
    
    tmp=[];
    cnt1=1;
    for j=1:size(b,2)
        tt=corr(b(i:i+11,j),b1(i:i+11,j));
        if ~isnan(tt)
            tmp(cnt1)=tt;
            cnt1=cnt1+1;
        end
    end
    if ~isempty(tmp)
        all_corr(2,cnt)=mean(tmp);
    end
    
    
    tmp=[];
    cnt1=1;
    for j=1:size(c,2)
        tt=corr(c(i:i+11,j),c1(i:i+11,j));
        if ~isnan(tt)
            tmp(cnt1)=tt;
            cnt1=cnt1+1;
        end
    end
    if ~isempty(tmp)
        all_corr(3,cnt)=mean(tmp);
    end
    cnt=cnt+1;
end
figure
plot(all_corr','linewidth',2)
legend('wt','cone','rod')
line([0,13],[0,0],'color','k')
set(gca,'XTick',1:12,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})



a=wt(:,:,2)';
b=cone(:,:,2)';
c=rod(:,:,2)';
a1=abs(wt(:,:,1)');
b1=abs(cone(:,:,1)');
c1=abs(rod(:,:,1)');

cnt=1;
figure
hold on
for i=1:12:144
    tmp=[];
    cnt1=1;
    for j=1:size(a,2)
        tt=corr(a(i:i+11,j),a1(i:i+11,j));
        if ~isnan(tt)
            tmp(cnt1)=tt;
            cnt1=cnt1+1;
        end
    end
    plot(ones(1,cnt1-1)*cnt,tmp,'b*')
    
    tmp=[];
    cnt1=1;
    for j=1:size(b,2)
        tt=corr(b(i:i+11,j),b1(i:i+11,j));
        if ~isnan(tt)
            tmp(cnt1)=tt;
            cnt1=cnt1+1;
        end
    end    
    plot(ones(1,cnt1-1)*cnt+0.2,tmp,'g*')
    
    
%     tmp=[];
%     cnt1=1;
%     for j=1:size(c,2)
%         tt=corr(c(i:i+11,j),c1(i:i+11,j));
%         if ~isnan(tt)
%             tmp(cnt1)=tt;
%             cnt1=cnt1+1;
%         end
%     end
%     plot(ones(1,cnt1-1)*cnt,tmp,'r*')
    cnt=cnt+1;
end
line([0,13],[0,0],'color','k')
set(gca,'XTick',1:12,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})



