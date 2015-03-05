clear
date='20120329'
date='20120627'
date='20120714'
date='20120928'
date='20121002'
date='20121023'
date='20130220'
date='20130220_1'
date='20130224'
date='20130225'
date='20130226'
date='20130227'
date='20130301'
date='20130301_1'
date='20130301_2'
date='20130302'
date='20130302_1'

date='20120902_1'
date='20120902_2'

date='20121023_1'
date='20121026_1'
codeWord='quick';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% select heka files
file_list=[];
ff=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
    end
    if ~isempty(regexp(heka(i).name,'nd_8', 'once'))
        ff=[ff i];
    end
end
if ~isempty(ff)
    file_list(file_list<ff(1))=[];
end

% get protocol times
protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol=round(protocol(2:end-2,[1,2]));
protocol(protocol(1:end,2)==1,2)=50;
protocol(protocol(1:end,2)==0,2)=10;
protocol(protocol(1:end,2)==-1,2)=30;

units=dir([mainpath,'units/*.mat']);
white_flash=zeros(4500,length(file_list),length(units));
black_flash=zeros(4500,length(file_list),length(units));
names=cell(length(units),1);
for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        conv=convolved(spikes,40,75000);
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash(:,i,cnt)=white_flash(:,i,cnt)+conv(protocol(st,1)-499:protocol(st,1)+4000)'/5;
            black_flash(:,i,cnt)=black_flash(:,i,cnt)+conv(protocol(st+2,1)-499:protocol(st+2,1)+4000)'/5;
        end
    end  
    
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
end

save([path2save,'quick'],'white_flash','black_flash','names')


%% Accumulate

clear

dates=cell(15,1);
dates{1}='20130220'
dates{2}='20130220_1'
dates{3}='20130224'
dates{4}='20130225'
dates{5}='20130226'
dates{6}='20130227'
dates{7}='20130301'
dates{8}='20130301_1'
dates{9}='20130301_2'
dates{10}='20130302'
dates{11}='20130302_1'
dates{12}='20120329'
dates{13}='20120627'
dates{14}='20120714'
dates{15}='20121023'

wfON=cell(8,1);
bfON=cell(8,1);
wfOFF=cell(8,1);
bfOFF=cell(8,1);
nameON=[];
nameOFF=[];

for datesCNT=1:15
    date=dates{datesCNT}
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
    load([path2save,'quick'],'white_flash','black_flash','names')
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'onOff')
    size(black_flash,2)
    nameON=[nameON; names(onOff>0)];
    nameOFF=[nameOFF; names(onOff<0)];
    t=1;
    if size(black_flash,2)==72
        for i=8:9:72   
            tmp=mean(white_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=mean(white_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==32
        for i=3:4:32
            tmp=mean(white_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=mean(white_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==16
        for i=2:2:16
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==17
        for i=2:2:16
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==24
        for i=3:3:24
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==56
        for i=4:5:40
            tmp=mean(white_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=mean(white_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==8
        for i=1:8
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==7
        t=2;
        for i=1:7
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
        wfON{1}=[wfON{1} zeros(4500,sum(onOff>0))];
        wfOFF{1}=[wfOFF{1} zeros(4500,sum(onOff<0))];
        bfON{1}=[bfON{1} zeros(4500,sum(onOff>0))];
        bfOFF{1}=[bfOFF{1} zeros(4500,sum(onOff<0))];
    elseif size(black_flash,2)==14
        for i=1:8
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    end    

end

clear black_flash white_flash date dates dateCNT tmp path2save names onOff




%% Plots

% Plot every cell separately
path2save='/mnt/muench_data/data/alexandra/MEA_data/analysis/Quick_OFF_cells/';
if ~exist(path2save,'dir')
    mkdir(path2save)
end
figure
set(gcf,'position',[20 50        1639         900])
nds='87654321'
for cnt=1:size(bfOFF{2},2)
    maxFR=[];
    for i=1:8
        subplot(4,2,i)
        hold off
        plot(bfOFF{i}(:,cnt),'r','LineWidth',2)
        hold on
        maxFR=[maxFR max(bfOFF{i}(:,cnt))];
        plot(4701:9200,wfOFF{i}(:,cnt),'LineWidth',2)
        maxFR=[maxFR max(wfOFF{i}(:,cnt))];        
    end
    k=max(maxFR)*1.1;
    for i=1:8
        subplot(4,2,i)
        axis([0 9200 0 k]);
        line([500,500],[0 k],'color','k')
        line([2500,2500],[0 k],'color','k')
        line([5200,5200],[0 k],'color','k')
        line([7200,7200],[0 k],'color','k')
        title(['ND',nds(i)])
    end
    subplot('Position',[0.5 0.96 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([nameON{cnt},'   RED - BLACK FLASH, BLUE - WHITE FLASH'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,nameON{cnt},'.png'])
end



a=bfON{4}(500:3500,:);
x=cov(a');
[V,~]=eig(x);
pc_vectors=V(:,end-2:end);
pc1=a'*pc_vectors(:,end);
pc2=a'*pc_vectors(:,end-1);
figure
plot(pc1,pc2,'.')
% plot mean of all cells (for ON and OFF), substract spontaneous FR for
% every cell

figure
nds='87654321';
for i=1:8
    subplot(4,2,i)
    hold on
    a=wfON{i};
%     a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(mean(a,2),'r','lineWidth',3)
    a=bfON{i};
%     a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(4701:9200,mean(a,2),'r','lineWidth',3)
    a=wfOFF{i};
%     a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(4701:9200,mean(a,2),'b','lineWidth',3)
    a=bfOFF{i};
%     a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(mean(a,2),'b','lineWidth',3)
    axis([0 9200 -15 55])
    line([500,500],[-15,55],'color','k')
    line([2500,2500],[-15,55],'color','k')
    line([5300,5300],[-15,55],'color','k')
    line([7300,7300],[-15,55],'color','k')
    line([0,9200],[0,0],'color','k')
    text(550,52,'Positive st')
    text(5350,52,'Negative st')
    title(['ND', nds(i), '  red ON, blue OFF'])
end






figure
cc=1;
for i=1:8
    subplot(8,2,cc)
    hold on
    a=wfON{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(a)
    plot(mean(a,2),'k','lineWidth',3)
    a=bfON{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(4501:9000,a)
    plot(4501:9000,mean(a,2),'k','lineWidth',3)
    
    subplot(8,2,cc+1)
    hold on
    a=wfOFF{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(a)
    plot(mean(a,2),'k','lineWidth',3)
    
    a=bfOFF{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(4501:9000,a)
    plot(4501:9000,mean(a,2),'k','lineWidth',3)
    
    cc=cc+2;
end



figure

i=4
hold off
a=bfOFF{i};
a=a-repmat(mean(a(1:490,:)),4500,1);
plot(a)
hold on
plot(mean(a,2),'k','lineWidth',3)
a=wfOFF{i};
a=a-repmat(mean(a(1:490,:)),4500,1);
plot(4701:9200,a)
plot(4701:9200,mean(a,2),'k','lineWidth',3)
axis([0 9200 -50 170])
line([500,500],[-50,170],'color','k','linewidth',2)
line([2500,2500],[-50,170],'color','k','linewidth',2)
line([5300,5300],[-50,170],'color','k','linewidth',2)
line([7300,7300],[-50,170],'color','k','linewidth',2)
line([0,9200],[0,0],'color','k','linewidth',2)
text(550,150,'Positive st')
text(5350,150,'Negative st')
title(['ND', nds(i), '  OFF cells n=',int2str(size(bfOFF{1},2))])


%for ON-OFF
line([2550 2950],[20,20],'color','r','linewidth',4)
line([550 950]+4700,[20,20],'color','r','linewidth',4)

line([2350 2450],[20,20],'color','c','linewidth',4)
line([350 450]+4700,[20,20],'color','c','linewidth',4)


%for rebound

line([3000 3400],[20,20],'color','r','linewidth',4)
line([1000 1400]+4700,[20,20],'color','r','linewidth',4)

line([2350 2450],[20,20],'color','c','linewidth',4)
line([350 450]+4700,[20,20],'color','c','linewidth',4)











% ON-OFF response
acc=[];
acc1=[];
for i=2:8
    a=bfOFF{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    b=mean(a(2550:2950,:));
    b1=mean(a(2350:2450,:));
    a=wfOFF{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    c=mean(a(550:950,:));
    c1=mean(a(350:450,:));
    d=b>b1&c>c1&b>0&c>0;
    acc1=[acc1 d'];


    a=wfON{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    b=mean(a(2550:2950,:));
    b1=mean(a(2350:2450,:));
    a=bfON{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    c=mean(a(550:950,:));
    c1=mean(a(350:450,:));
    d=b>b1&c>c1&b>0&c>0;
    acc=[acc d'];
end


figure
acc=logical(acc)
for i=2:8
subplot(2,4,i)
    hold on
    a=wfON{i}(:,acc(:,i-1));
    sum(acc(:,i-1))
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(a)
    plot(mean(a,2),'k','lineWidth',3)
    a=bfON{i}(:,acc(:,i-1));
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(4701:9200,a)
    plot(4701:9200,mean(a,2),'k','lineWidth',3)
    title([num2str(sum(acc(:,i-1))),'    ', num2str(sum(acc(:,i-1))/size(acc,1))])
    line([500,500],[-25,255],'color','k')
    line([2500,2500],[-25,255],'color','k')
    line([5300,5300],[-25,255],'color','k')
    line([7300,7300],[-25,255],'color','k')
    line([0,9200],[0,0],'color','k')
    axis([0 9200 -25 150])
end


figure
acc1=logical(acc1)
for i=2:8
subplot(2,4,i)
    hold on
    a=bfOFF{i}(:,acc1(:,i-1));
    sum(acc1(:,i-1))
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(a)
    plot(mean(a,2),'k','lineWidth',3)
    a=wfOFF{i}(:,acc1(:,i-1));
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(4701:9200,a)
    plot(4701:9200,mean(a,2),'k','lineWidth',3)
    title([num2str(sum(acc1(:,i-1))),'    ', num2str(sum(acc1(:,i-1))/size(acc1,1))])
    line([500,500],[-25,255],'color','k')
    line([2500,2500],[-25,255],'color','k')
    line([5300,5300],[-25,255],'color','k')
    line([7300,7300],[-25,255],'color','k')
    line([0,9200],[0,0],'color','k')
    axis([0 9200 -25 150])
end

figure
bar([sum(acc1)/size(acc1,1)*100; sum(acc)/size(acc,1)*100]')
legend({'OFF','ON'})
set(gca,'xtick',1:7,'xticklabel',{'7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('% of total')
title('portion of cells with ON-OFF response at each ND')

figure
[a,~]=hist([sum(acc1,2);9]);
[b,~]=hist([sum(acc,2);9]);
plot(0:7,a(1:8)/size(acc,1)*100,'r')
bar([a(1:8)/size(acc1,1)*100; b(1:8)/size(acc,1)*100]')
legend('OFF','ON')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('% of total')
title('portion of cells with ON-OFF response at each ND')

% rebound response
acc=[];
acc1=[];
for i=2:8
    a=bfOFF{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    b=mean(a(3000:3400,:));
    b1=mean(a(2350:2450,:));
    a=wfOFF{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    c=mean(a(1000:1400,:));
    c1=mean(a(350:450,:));
    d=b>b1&c>c1&b>0&c>0;
    acc1=[acc1 d'];


    a=wfON{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    b=mean(a(3000:3400,:));
    b1=mean(a(2350:2450,:));
    a=bfON{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    c=mean(a(1000:1400,:));
    c1=mean(a(350:450,:));
    d=b>b1&c>c1&b>0&c>0;
    acc=[acc d'];
end


figure
acc=logical(acc)
for i=2:8
subplot(2,4,i)
    hold on
    a=wfON{i}(:,acc(:,i-1));
    sum(acc(:,i-1))
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(a)
    plot(mean(a,2),'k','lineWidth',3)
    a=bfON{i}(:,acc(:,i-1));
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(4701:9200,a)
    plot(4701:9200,mean(a,2),'k','lineWidth',3)
    title([num2str(sum(acc(:,i-1))),'    ', num2str(sum(acc(:,i-1))/size(acc,1))])
    line([500,500],[-25,255],'color','k')
    line([2500,2500],[-25,255],'color','k')
    line([5300,5300],[-25,255],'color','k')
    line([7300,7300],[-25,255],'color','k')
    line([0,9200],[0,0],'color','k')
    axis([0 9200 -25 150])
end


figure
acc1=logical(acc1)
for i=2:8
subplot(2,4,i)
    hold on
    a=bfOFF{i}(:,acc1(:,i-1));
    sum(acc1(:,i-1))
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(a)
    plot(mean(a,2),'k','lineWidth',3)
    a=wfOFF{i}(:,acc1(:,i-1));
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(4701:9200,a)
    plot(4701:9200,mean(a,2),'k','lineWidth',3)
    title([num2str(sum(acc1(:,i-1))),'    ', num2str(sum(acc1(:,i-1))/size(acc1,1))])
    line([500,500],[-25,255],'color','k')
    line([2500,2500],[-25,255],'color','k')
    line([5300,5300],[-25,255],'color','k')
    line([7300,7300],[-25,255],'color','k')
    line([0,9200],[0,0],'color','k')
    axis([0 9200 -25 150])
end



figure
bar([sum(acc1)/size(acc1,1)*100; sum(acc)/size(acc,1)*100]')
legend({'OFF','ON'})
set(gca,'xtick',1:7,'xticklabel',{'7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('% of total')
title('portion of cells with rebound at each ND')

figure
[a,b]=hist([sum(acc1(:,1:7),2);9]);
plot(0:7,a(1:8)/size(acc1,1)*100)
hold on
[a,b]=hist([sum(acc(:,1:7),2);9]);
plot(0:7,a(1:8)/size(acc,1)*100,'r')

sum(sum(acc1(:,[3 5]),2)==2)










figure
nds='87654321';
for i=1:8
    subplot(4,2,i)
    hold on
    a=wfON{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(mean(a,2),'r','lineWidth',3)
    a=bfON{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(500:2500,mean(a(2500:4500,:),2),'b','lineWidth',3)
    plot(2501:4500,mean(a(501:2500,:),2),'b','lineWidth',3)

    axis([0 4500 -15 50])
    title(['ND', nds(i), '  red white flash, blue black flash, spont substr'])
end


figure
nds='87654321';
for i=1:8
    subplot(4,2,i)
    hold on
    a=wfOFF{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(mean(a,2),'r','lineWidth',3)
    a=bfOFF{i};
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(500:2500,mean(a(2500:4500,:),2),'b','lineWidth',3)
    plot(2501:4500,mean(a(501:2500,:),2),'b','lineWidth',3)

    axis([0 4500 -15 50])
    title(['ND', nds(i), ' OFF'])
end





col='brgmc';
for i=1:size(bfOFF{2},2)
    if mod(i-1,49)==0
        figure
        cc=1;
    end
    subplot(7,7,cc)
    hold on
    for j=2:6
        a=bfOFF{j}(:,i);
%         a=a-mean(a(1:490));
        plot(a,col(j-1))
    end
    set(gca,'Xtick',0,'xticklabel','')
    set(gca,'ytick',0,'yticklabel','')
    line([0,4500],[0,0],'color','k')
    axis tight
    a=get(gca,'Ylim');
    line([2500,2500],[a(1),a(2)],'color','k')
    line([500,500],[a(1),a(2)],'color','k')
    title(int2str(i))

    cc=cc+1;
end





for j=2:8
    for i=1:size(bfOFF{2},2)
        a=bfOFF{j}(:,i);
        b=wfOFF{j}(:,i);
        a=a-mean(a(1:490));
%         transOFF(i,j)=sum(a(500:850))/sum(b(2500:2850));
%         transOFF(i,j)=sum(abs(a(500:850)))/sum(abs(a(2150:2500)));
        transOFF(i,j)=sum(a(500:850))-sum(a(2900:3500));
    end
    for i=1:size(wfON{2},2)
        a=wfON{j}(:,i);
        b=bfON{j}(:,i);
        a=a-mean(a(1:490));
%         transON(i,j)=sum(a(500:850))/sum(b(2500:2850));
%         transON(i,j)=sum(abs(a(500:850)))/sum(abs(a(2150:2500)));
        transON(i,j)=sum(a(500:850))-sum(a(2900:3500));
    end
end


figure
plot(transOFF(:,3),transOFF(:,5),'.b')
hold on
plot(transON(:,3),transON(:,5),'.r')
line([-100 100]*1000,[-100 100]*1000,'color','k')
figure
subplot(1,2,1)
plot(transOFF')
subplot(1,2,2)
plot(transON')
figure
hist(transOFF(:,3),30)


figure;
clear a b c
m=wfOFF;
for i=1:size(m{2},2)
    a(i)=corr(m{3}(500:2500,i),m{4}(500:2500,i));
    b(i)=corr(m{3}(500:2500,i),m{5}(500:2500,i));
    c(i)=corr(m{4}(500:2500,i),m{5}(500:2500,i));
end
[~, k]=sort(b);
plot([a(k); b(k); c(k)]')


