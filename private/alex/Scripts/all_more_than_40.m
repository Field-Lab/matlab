


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


subplot(3,1,2)
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[-400, 400],'color','k')    
end


a1=abs(wt(:,:,1)');
b1=abs(cone(:,:,1)');
c1=abs(rod(:,:,1)');
clear k
for i=1:size(a1,1)
    tmp=a1(i,:)>40;
    k(i,1)=mean(a1(i,tmp));
    tmp=b1(i,:)>40;
    k(i,2)=mean(b1(i,tmp));
    tmp=c1(i,:)>40;
    k(i,3)=mean(c1(i,tmp));
end
k(isnan(k))=0;
plot(k,'.','MarkerSize',10)

line([0,147],[0,0],'color','k')
axis([0,147,-400,400])
title('Amplitude')


subplot(3,1,1)
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})

a=wt(:,:,2)';
b=cone(:,:,2)';
c=rod(:,:,2)';
clear k
for i=1:size(a,1)
    tmp=a(i,:)>0&a1(i,:)>40;
    k(i,1)=mean(a(i,tmp));
    tmp=b(i,:)>0&b1(i,:)>40;
    k(i,2)=mean(b(i,tmp));
    tmp=c(i,:)>0&c1(i,:)>40;
    k(i,3)=mean(c(i,tmp));
end
k(isnan(k))=0;

plot(k,'.','MarkerSize',10)
legend('wt', 'cone','rod')

for i=12.5:12:146
    line([i,i],[0, 200],'color','k')    
end
axis([0,147,50,200])
title('Latency')



subplot(3,1,3)
hold on
set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
for i=12.5:12:146
    line([i,i],[0, 80],'color','k')    
end

a=wt(:,:,3)';
b=cone(:,:,3)';
c=rod(:,:,3)';
clear k
for i=1:size(a,1)
for i=1:size(a,1)
    tmp=a(i,:)>0&a1(i,:)>40;
    k(i,1)=mean(a(i,tmp));
    tmp=b(i,:)>0&b1(i,:)>40;
    k(i,2)=mean(b(i,tmp));
    tmp=c(i,:)>0&c1(i,:)>40;
    k(i,3)=mean(c(i,tmp));
end
end
k(isnan(k))=0;

plot(k,'.','MarkerSize',10)
line([0,147],[0,0],'color','k')
axis([0,147,0,80])
title('Width')
























figure
set(gcf,'position',[1 31 1680 946])

subplot(3,1,2)
hold on
set(gca,'XTick',6:12:36,'xticklabel',{'4','3','2'})
for i=12.5:12:36
    line([i,i],[-400, 400],'color','k')    
end


a1=abs(wt(:,49:84,1)');
b1=abs(cone(:,49:84,1)');
c1=abs(rod(:,49:84,1)');
clear k
for i=1:size(a1,1)
    tmp=a1(i,:)>40;
    k(i,1)=mean(a1(i,tmp));
    tmp=b1(i,:)>40;
    k(i,2)=mean(b1(i,tmp));
    tmp=c1(i,:)>40;
    k(i,3)=mean(c1(i,tmp));
end
k(isnan(k))=0;
plot(k,'.','MarkerSize',10)

line([0,37],[0,0],'color','k')
axis([0,37,-400,400])
title('Amplitude')


subplot(3,1,1)
hold on
set(gca,'XTick',6:12:36,'xticklabel',{'4','3','2'})

a=wt(:,49:84,2)';
b=cone(:,49:84,2)';
c=rod(:,49:84,2)';
clear k
for i=1:size(a,1)
    tmp=a(i,:)>0&a1(i,:)>40;
    k(i,1)=mean(a(i,tmp));
    tmp=b(i,:)>0&b1(i,:)>40;
    k(i,2)=mean(b(i,tmp));
    tmp=c(i,:)>0&c1(i,:)>40;
    k(i,3)=mean(c(i,tmp));
end
k(isnan(k))=0;

plot(k,'.','MarkerSize',10)
legend('wt', 'cone','rod')

for i=12.5:12:36
    line([i,i],[0, 200],'color','k')    
end
axis([0,37,50,200])
title('Latency')



subplot(3,1,3)
hold on
set(gca,'XTick',6:12:36,'xticklabel',{'4','3','2'})
for i=12.5:12:36
    line([i,i],[0, 80],'color','k')    
end

a=wt(:,49:84,3)';
b=cone(:,49:84,3)';
c=rod(:,49:84,3)';
clear k
for i=1:size(a,1)
for i=1:size(a,1)
    tmp=a(i,:)>0&a1(i,:)>40;
    k(i,1)=mean(a(i,tmp));
    tmp=b(i,:)>0&b1(i,:)>40;
    k(i,2)=mean(b(i,tmp));
    tmp=c(i,:)>0&c1(i,:)>40;
    k(i,3)=mean(c(i,tmp));
end
end
k(isnan(k))=0;

plot(k,'.','MarkerSize',10)
line([0,37],[0,0],'color','k')
axis([0,37,0,80])
title('Width')





