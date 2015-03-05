
figure
kk=1;
i=5;
cnt=1;
while kk<=15

    if onOff_103(cnt)<0
%         b=(round(correlFiringRate(i,cnt)*100))/100;
%         if b<-0.2
            subplot(5,3,kk)
            if onOff_103(cnt)>0
                plot(FiringRate(1:600000,i,cnt),'r')
            else
                plot(FiringRate(1:600000,i,cnt),'b')
            end
            title([int2str(cnt),'   ',int2str(onOff(cnt))])
            kk=kk+1;
%         end
    end
    cnt=cnt+1;
end


figure
kk=1;
i=5;
cnt=33;
while kk<=15

    if onOff_103(cnt)>0
%         b=(round(correlFiringRate(i,cnt)*100))/100;
%         if b<-0.2
            subplot(5,3,kk)
            if onOff_103(cnt)>0
                plot(FiringRate(1:600000,i,cnt),'r')
            else
                plot(FiringRate(1:600000,i,cnt),'b')
            end
            title([int2str(cnt),'   ',int2str(onOff(cnt))])
            kk=kk+1;
%         end
    end
    cnt=cnt+1;
end




figure
kk=1;
i=5;
cnt=65;
while kk<=15

    if onOff_103(cnt)>0
%         b=(round(correlFiringRate(i,cnt)*100))/100;
%         if b<-0.2
            subplot(5,3,kk)
            if onOff_103(cnt)>0
                plot(FiringRate(1:600000,i,cnt),'r')
            else
                plot(FiringRate(1:600000,i,cnt),'b')
            end
            title([int2str(cnt),'   ',int2str(onOff(cnt))])
            kk=kk+1;
%         end
    end
    cnt=cnt+1;
end





figure
kk=1;
i=5;
cnt=1;
while kk<=25
    
    if onOff_103(cnt)>0
        subplot(5,5,kk)
        p=1;
        for j=1:500:600000-500
            ss(p)=std(FiringRate(j:j+1000,i,cnt));
            p=p+1;
        end
        plot(ss,'r')
        axis([0 1200 0 30])
        title([int2str(cnt),'   ',int2str(onOff(cnt))])
        kk=kk+1;
    end
    cnt=cnt+1;
end












%% 20130301
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130301'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
isi=zeros(200,size(SpikeCount,1),6,size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            for mc=1:6
                take=correctedProtocols(600*(mc-1)+1,1,i)+2000:correctedProtocols(600*mc,1,i)-500;
                spikesTMP=spikes(spikes>(correctedProtocols(600*(mc-1)+1,1,i)+1000)&spikes<(correctedProtocols(600*mc,1,i)-500));
                a=diff(spikesTMP);
                a(a>200)=[];
                if length(a)>4
                    a(a==0)=1;
                    [z,x]=hist(a,unique(a));
                    isi(x,i,mc,cnt)=isi(x,i,mc,cnt)+z';
                end
            end
        end
    end
end

figure
set(gcf,'position',[560          17        1677         931])
nds='87654321'
coord=sort([1:4:4*8 2:4:4*8]);
cc=1;
clear axislimits
for cnt=1:size(FiringRate,3)
    if onOff(cnt)>0
        col='r';
        col2='m';
    else
        col='b';
        col2='c';
    end
    if onOff(cnt)>0
        subplot(6,3,cc)
        for i=1:16
            tmp=FiringRate(:,i,cnt);
            for mc=1:6
                take=correctedProtocols(600*(mc-1)+1,1,i)+2000:correctedProtocols(600*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:7000));
            end
            tmp=reshape(tmp1(:,1:2:end),21000,1);
            ss=std(tmp);
            mm=mean(tmp);
            errorbar(coord(i),mm,ss,col)
            hold on
            plot(coord(i),mm,'.','markersize',20,'color',col)
            tmp=reshape(tmp1(:,2:2:end),21000,1);
            ss=std(tmp);
            mm=mean(tmp);
            errorbar(coord(i)+0.3,mm,ss,col2)
            plot(coord(i)+0.3,mm,'.','markersize',20,'color',col2)
        end
        line([0,31],[0,0],'color','k')
        axis tight
        axislimits(cc,:)=get(gca,'Ylim');
        for i=1:16
            line([3.5+(i-1)*4,3.5+(i-1)*4],[-100 100],'color','k')
        end
        axis([0 31 axislimits(cc,1)*1.1 axislimits(cc,2)*1.1])
        set(gca,'xtick',2:4:30,'xticklabel',{'8','7','6','5','4','3','2','1'})
        cc=cc+1;
    end
end

for cc=1:18
    subplot(6,3,cc)
    axis([0 31 min(axislimits(:,1)) max(axislimits(:,2))])
    hold off
end




%% 20121023
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121023'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
        end
    end
end

figure
set(gcf,'position',[560          17        1677         931])
nds='87654321'
coord=sort([1:4:4*8 2:4:4*8]);
cc=1;
clear axislimits
for cnt=1:size(FiringRate,3)
    if onOff(cnt)>0
        col='r';
        col2='m';
    else
        col='b';
        col2='c';
    end
    if onOff(cnt)>0
        subplot(6,5,cc)
        km=1;
        for i=1:24:24*8
            tmp=FiringRate(:,i:2:(i+23),cnt);
            tmp=reshape(tmp(2501:62500,:),720000,1);
            ss=std(tmp);
            mm=mean(tmp);
            errorbar(km,mm,ss,col)
            hold on
            plot(km,mm,'.','markersize',20,'color',col)
            tmp=FiringRate(:,(i+1):2:(i+23),cnt);
            tmp=reshape(tmp(2501:62500,:),720000,1);
            ss=std(tmp);
            mm=mean(tmp);
            errorbar(km+0.2,mm,ss,col2)
            plot(km+0.2,mm,'.','markersize',20,'color',col2)
            km=km+1;
        end
        line([0.5,8.7],[0,0],'color','k')
        axis tight
        axislimits(cc,:)=get(gca,'Ylim');
        for i=1:8
            line([1.7+(i-1),1.7+(i-1)],[-100 100],'color','k')
        end
        axis([0.5 8.7 axislimits(cc,1)*1.1 axislimits(cc,2)*1.1])
        set(gca,'xtick',1.2:8.2,'xticklabel',{'8','7','6','5','4','3','2','1'})
        title(names{cnt})
        cc=cc+1;
    end
end

for cc=1:29
    subplot(6,5,cc)
    axis([0.5 8.7 min(axislimits(:,1)) max(axislimits(:,2))])
    hold off
end



%% 20121002
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120928'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
        end
    end
end

figure
set(gcf,'position',[560          17        1677         931])
nds='76543212'
coord=sort([1:4:4*8 2:4:4*8]);
cc=1;
clear axislimits
for cnt=1:size(FiringRate,3)
    if onOff(cnt)>0
        col='r';
        col2='m';
    else
        col='b';
        col2='c';
    end
    if onOff(cnt)>0
        subplot(7,3,cc)
        km=1;
        for i=1:24:24*8
            tmp=FiringRate(:,i:2:(i+23),cnt);
            tmp=reshape(tmp(2501:62500,:),720000,1);
            ss=std(tmp);
            mm=mean(tmp);
            errorbar(km,mm,ss,col)
            hold on
            plot(km,mm,'.','markersize',20,'color',col)
            tmp=FiringRate(:,(i+1):2:(i+23),cnt);
            tmp=reshape(tmp(2501:62500,:),720000,1);
            ss=std(tmp);
            mm=mean(tmp);
            errorbar(km+0.2,mm,ss,col2)
            plot(km+0.2,mm,'.','markersize',20,'color',col2)
            km=km+1;
        end
        line([0.5,8.7],[0,0],'color','k')
        axis tight
        axislimits(cc,:)=get(gca,'Ylim');
        for i=1:8
            line([1.7+(i-1),1.7+(i-1)],[-100 100],'color','k')
        end
        axis([0.5 8.7 axislimits(cc,1)*1.1 axislimits(cc,2)*1.1])
        set(gca,'xtick',1.2:8.2,'xticklabel',{'7','6','5','4','3','2','1','2'})
        cc=cc+1;
    end
end

for cc=1:21
    subplot(7,3,cc)
    axis([0.5 8.7 min(axislimits(:,1)) max(axislimits(:,2))])
    hold off
end


%% SOMETHING ELSE



tmp=FiringRate(:,i,cnt);
for mc=1:6
    take=correctedProtocols(600*(mc-1)+1,1,i)+2000:correctedProtocols(600*mc,1,i)-500;
    tmp1(:,mc)=tmp(take(1:7000));
end
tmp=reshape(tmp1(:,1:2:end),21000,1);
subplot(2,1,1)
hist([tmp; 0; 40],200)
axis([0 40 0 2000])

tmp=reshape(tmp1(:,2:2:end),21000,1);
subplot(2,1,2)
hist([tmp; 0; 40],200)
axis([0 40 0 1500])


cnt=4
i=5
figure
tmp=FiringRate(1:correctedProtocols(1,1,i),i,cnt);
subplot(3,1,1)
hist([tmp; 0; 40],100)
% axis([0 40 0 700])

tmp=FiringRate(:,i,cnt);
for mc=1:6
    take=correctedProtocols(600*(mc-1)+1,1,i)+2000:correctedProtocols(600*mc,1,i)-500;
    tmp1(:,mc)=tmp(take(1:7000));
end
tmp=reshape(tmp1(:,1:2:end),21000,1);
subplot(3,1,2)
hist([tmp; 0; 40],100)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w')
axis([0 40 0 700])
hold on
tmp=FiringRate(1:correctedProtocols(1,1,i),i,cnt);
hist([tmp; 0; 40],100)




tmp=reshape(tmp1(:,2:2:end),21000,1);
subplot(3,1,3)
hist([tmp; 0; 40],100)
axis([0 40 0 700])


cnt=3
kk=1;
for i=5:2:12    
    subplot(2,4,kk)
    plot(sum(isi(1:150,i,1:2:end,cnt),3))
    subplot(2,4,kk+4)
    plot(sum(isi(1:150,i,2:2:end,cnt),3))
    kk=kk+1;
end