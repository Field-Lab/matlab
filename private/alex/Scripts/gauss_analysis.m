clear
path2coef='S:\user\alexandra\MEA_data\20120925\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,1:149,1:3)=common_res_fit';
end

subplot(2,1,1)
plot(res(:,13:120,2)','b.')
set(gca,'XTick',6:12:12*9,'xticklabel',{'5','4','3','2','3','2','3','2','1','2'})
for i=12.5:12:12*9
    line([i,i],[0, 200],'color','m')    
end
a=res(:,13:120,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'k.','MarkerSize',20)
axis([0,109,0,200])
subplot(2,1,2)
plot(abs(res(:,13:120,1))','b.')
set(gca,'XTick',6:12:12*9,'xticklabel',{'5','4','3','2','3','2','3','2','1','2'})
for i=12.5:12:12*9
    line([i,i],[0, 200],'color','m')    
end
a=abs(res(:,13:120,1))';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(abs(a(i,tmp)));
end
hold on
plot(b,'k.','MarkerSize',20)
axis([0,109,0,200])

a=res(:,13:120,1)';
b=res(:,13:120,2)';
cnt=1;
nd='543232321';
for i=1:12:12*9
    subplot(3,3,cnt)
    t=diff(mean(abs(a(i:i+11,:))));
    tt=diff(mean(b(i:i+11,:)));
    plot(t,tt,'.')
    title(['ND',nd(cnt)])
    cnt=cnt+1;
end
plot(a,b,'.')







clear
path2coef='S:\user\alexandra\MEA_data\20121004\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,:,1:3)=common_res_fit';
end
figure
subplot(3,1,1)
plot(res(:,:,2)','b.')
% set(gca,'XTick',6:12:12*14,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
% for i=12.5:12:12*14
%     line([i,i],[0, 200],'color','m')    
% end
a=res(:,:,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'k.','MarkerSize',20)
axis([0,173,0,200])
title('Latency')

subplot(3,1,2)
plot(res(:,:,1)','b.')
set(gca,'XTick',6:12:12*14,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*14
    line([i,i],[-700, 700],'color','m')    
end
line([0,121],[0,0],'color','k')
axis([0,173,-700,700])
title('Amplitude')

subplot(3,1,3)
plot(res(:,:,3)','b.')
set(gca,'XTick',6:12:12*14,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*14
    line([i,i],[0, 80],'color','m')    
end
line([0,173],[0,0],'color','k')
axis([0,173,0,80])
title('Width')





figure
subplot(3,1,1)
bar(mean(res(:,1:172,2)))
set(gca,'XTick',6:12:12*14,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*14
    line([i,i],[0, 150],'color','m')    
end
hold on
title('Latency')

subplot(3,1,2)
bar(mean(abs(res(:,1:172,1))))
set(gca,'XTick',6:12:12*14,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*14
    line([i,i],[0, 350],'color','m')    
end
title('Amplitude')

subplot(3,1,3)
bar(mean(res(:,1:172,3)))
set(gca,'XTick',6:12:12*14,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*14
    line([i,i],[0, 50],'color','m')    
end
title('Width')







clear
path2coef='S:\user\alexandra\MEA_data\20120921\fit_res_HC\';
coef=dir([path2coef,'*.mat']);
load([path2coef,coef(1).name])
toTr=length(common_res_fit);
for i=1:length(coef)
    load([path2coef,coef(i).name])
    res(i,1:toTr,1:3)=common_res_fit';
end

subplot(3,1,1)
plot(res(:,1:toTr,2)','b.')
set(gca,'XTick',6:12:12*(toTr/12+1),'xticklabel',{'7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*(toTr/12+1)
    line([i,i],[0, 200],'color','m')    
end
a=res(:,1:toTr,2)';
for i=1:size(a,1)
    tmp=find(a(i,:));
    b(i)=mean(a(i,tmp));
end
hold on
plot(b,'k.','MarkerSize',20)
axis([0,toTr,0,200])
title('Latency')

subplot(3,1,2)
plot(res(:,1:toTr,1)','b.')
set(gca,'XTick',6:12:12*(toTr/12+1),'xticklabel',{'7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*(toTr/12+1)
    line([i,i],[-700, 700],'color','m')    
end
line([0,toTr],[0,0],'color','k')
axis([0,toTr,-700,700])
title('Amplitude')

subplot(3,1,3)
plot(res(:,1:toTr,3)','b.')
set(gca,'XTick',6:12:12*(toTr/12+1),'xticklabel',{'7','6','5','4','3','2','1','2','3','4','5','6','7'})
for i=12.5:12:12*(toTr/12+1)
    line([i,i],[0, 80],'color','m')    
end
line([0,toTr],[0,0],'color','k')
axis([0,toTr,0,80])
title('Width')


%% Furier
load('F:\20121023\FFFlicker_LF\A20121023_CH24_sort1_100_unit_0009_FFFlicker_linear_filter.mat')
b=zeros(172,257);
for i=1:172
    filterEx=HighFilters(:,i);
    
    m=filterEx;
    NFFT = 2^nextpow2(numel(m));
    Y = fft(m,NFFT)/numel(m);
    f = 10000/2*linspace(0,1,NFFT/2+1);
    b(i,:)=2*abs(Y(1:NFFT/2+1));
end
figure
for i=1:172
    val=cumsum(b(i,2:20));
    val=val(end)/3;
    a(i)=find(cumsum(b(i,2:20))>=val,1);
end
plot(a)
for i=12.5:12:12*14
    line([i,i],[0, 20],'color','m')    
end
figure
for i=1:12*8
    subplot(8,12,i)
    plot(f(1:20),b(i,1:20))
end
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
