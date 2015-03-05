%% Fig.O
clear

date='20130302_1'
codeWord='NatMov'
dataRaster=cell(3,1);
dataRaster_supp=zeros(3,1);
data=zeros(4500,2,3);

mainpath=['/Users/alexth/Desktop/old_stuff/my_data/',date,'/'];
hekapath=['/Users/alexth/Desktop/old_stuff/my_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);


% select heka files
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
    end
end
conv_acc=zeros(25000,18,3);
units=dir([mainpath, 'units/*.mat']);
cnt=26;
ordCnt=387;

cnt=27
units(cnt).name
load([mainpath,'units/',units(cnt).name]);

t=[];

cnt=0;
for i=37:37+17
    spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms  
    
%     tmp=rand(100,1);
%     tmp=round(tmp*25000);
%     spikes=sort([spikes tmp']);
    conv_tmp=convolved(spikes,40,25000);
    conv_acc(:,cnt+1,1)=conv_tmp(121:end-120)/18;
    t=[t spikes+26000*cnt];    
    cnt=cnt+1;
end

dataRaster{1}=t;
dataRaster_supp(1)=cnt;

t=[];
cnt=0;
for i=55:55+17
    spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms 

    conv_tmp=convolved(spikes,40,25000);
    conv_acc(:,cnt+1,2)=conv_tmp(121:end-120)/18;
    t=[t spikes+26000*cnt];    
    cnt=cnt+1;
end

dataRaster{2}=t;
dataRaster_supp(2)=cnt;

t=[];
cnt=0;
for i=73:73+17
    spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
    conv_tmp=convolved(spikes,40,25000);
    conv_acc(:,cnt+1,3)=conv_tmp(121:end-120)/18;
    t=[t spikes+26000*cnt];    
    cnt=cnt+1;
end

dataRaster{3}=t;
dataRaster_supp(3)=cnt;

figure
nds='654';
cnt=1;
for i=1:3
    h=subplot(3,2,cnt)
    rasterplot(dataRaster{1},dataRaster_supp(2),26000,h);
    title(['ND',nds(i)])
    subplot(3,2,cnt+1)
    plot(mean(conv_acc(:,:,i),2))
    title(['ND',nds(i)])
    cnt=cnt+2;
    axis tight
end




steps=200;
ranges=1:100:25000-steps;
allCorr=zeros(numel(ranges),3);
cnt=1;
for j=ranges
    myCorr=triu(corr(conv_acc(j:j+steps,:,1)),1);
    myCorr(myCorr==0)=NaN;
    tmp=nanmean(myCorr);
    if sum(isnan(tmp))<size(conv_acc,2)/2% && nanmean(tmp)>0.7
        allCorr(cnt,1)=nanmean(tmp);        
    end
    
    myCorr=triu(corr(conv_acc(j:j+steps,:,2)),1);
    myCorr(myCorr==0)=NaN;
    tmp=nanmean(myCorr);
    if sum(isnan(tmp))<size(conv_acc,2)/2 %&& nanmean(tmp)>0.7
        allCorr(cnt,2)=nanmean(tmp);        
    end
    
    myCorr=triu(corr(conv_acc(j:j+steps,:,3)),1);
    myCorr(myCorr==0)=NaN;
    tmp=nanmean(myCorr);
    if sum(isnan(tmp))<size(conv_acc,2)/2 %&& nanmean(tmp)>0.7
        allCorr(cnt,3)=nanmean(tmp);        
    end
    
    cnt=cnt+1;
end




figure
plot(ranges,allCorr)

axis tight
hold on
plot(ranges,allCorr(:,1:2)*25)

tmp=allCorr(:,1)-allCorr(:,2);
plot(tmp)





% according to Thomas
myCorr=triu(corr(a),1);
myCorr(myCorr==0)=NaN;
myCorr=myCorr(:);
myCorr(isnan(myCorr))=[];
[x,y]=hist(myCorr,25);


myCorr=triu(corr(b),1);
myCorr(myCorr==0)=NaN;
myCorr=myCorr(:);
myCorr(isnan(myCorr))=[];
[x2,y2]=hist(myCorr,25);


for i=1:18
    for j=1:18
        m(i,j)=corr(a(:,i),b(:,j));
    end
end
[x1,y1]=hist(m(:),25);

figure
plot(y,x/153)
hold on
plot(y1,x1/324,'r')
plot(y2,x2/153,'g')

figure
line([0.92-0.052,0.92+0.052],[1,1])
line([0.89-0.068,0.89+0.068],[1.25,1.25],'color','r')
line([0.91-0.047,0.91+0.047],[1.5,1.5],'color','g')
axis([0.5 1 0.5 2])


figure
line([0.7-0.079,0.7+0.079],[1,1])
line([0.62-0.12,0.62+0.12],[1.25,1.25],'color','r')
line([0.52-0.12,0.52+0.12],[1.5,1.5],'color','g')
axis([0.5 1 0.5 2])


figure

plot(diff(spikes))
hold on
plot(diff(spikes1),'r')


m=zeros(27000,18);

cnt=1;
for i=73:73+17
    spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
    m(spikes,cnt)=1;
    cnt=cnt+1;
end
k=[];
for i=1:size(m,1)-5
    k(i)=nnz(m(i:i+5,:))/18;
end
    plot(k)
m5=squeeze(sum(reshape(m,5400,5,18),2));

plot(sum(m5'))

