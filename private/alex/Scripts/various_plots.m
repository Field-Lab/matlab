%% FILTERS
% HIGH PASS - for SPIKES
passH=500;%frequency of the high-pass filter
sf=25000;%sampling frequency
fNormH = passH/(sf/2);
orderH=12;%filter order
[aH,bH] = butter(orderH, fNormH, 'high');
   
% LOW PASS - for LFP - 180 Hz
sf=25000;%sampling frequency
passL=180;orderL=2;
fNormL = passL/(sf/2);
[aL,bL] = butter(orderL, fNormL, 'low');

   
% LOW PASS - for LFP - 30 Hz
sf=25000;%sampling frequency
passL=30;orderL=10;
fNormL = passL/(sf/2);
[aL,bL] = butter(orderL, fNormL, 'low');


%% PLOTS
% plot lfp by channel from 0805
for ch=53:57
    if ch~=15
        load(['F:\20110818\processing\preprocessed\lfp\quick_CH',int2str(ch),'_0818']);
        close all
        figure(1)
        set(gcf,'Position',[1  31 1680 946]);
        colors=[0 0 255; 255 0 0; 0 255 0; 255 0 255; 200 170 50; 0 255 255; 120 50 200; 0 0 0]/255;
        sf=250;%sampling frequency
        passL=30;orderL=10;
        fNormL = passL/(sf/2);
        [aL,bL] = butter(orderL, fNormL, 'low');
        numbers=[2         151         300         449         598         747         896        1045        1194];
        leg='1234567';
        abs_min=0; abs_max=0;
        xlab=2.4:0.004:7.2;
        for j=1:7
            cnt=numbers(j);
            for i=0:20:130
                subplot(3,3,j)
                tmp=filtfilt(aL, bL, mean(lfp(6000:10:18000,cnt+i:2:cnt+i+19),2));
                plot(xlab,tmp,'color',colors(i/20+1,:))
                hold on
                abs_min=min(abs_min,min(min(tmp)));
                abs_max=max(abs_max,max(max(tmp)));
            end
            title(['ND ',int2str(10-j)])
            xlabel('s')
            ylabel('a.u.')
            if j==1
                legend(leg')
            end
        end
        for j=1:7
            subplot(3,3,j)
            axis([2.4 7.2 abs_min abs_max])
        end
        % subplot('Position',[0.2 0.97 0.6 0.02])
        % set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        % text(0.02,0.5,['20110818  CH',int2str(ch),', lfp <30Hz, ND 9-1, '])
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(['20110818  CH',int2str(ch),', LFP, low pass filt. 30Hz, boxes: ND 9-1, lines in boxes: mean of each 10 subs. trials of CONTRAST 2'],'FontSize',15,'FontWeight','bold')
        saveas(gcf,['F:\20110818\processing\pictures\lfp_by_channel\20110818_CH',int2str(ch),'.bmp']);
    end
end

% plot spikes by channel minus mean spont from 0805
for ch=1:60
    if ch~=15
        load(['F:\20110818\processing\preprocessed\spikes\quick_CH',int2str(ch),'_0818']);
        close all
        figure(1)
        set(gcf,'Position',[1  31 1680 946]);
        colors=[0 0 255; 255 0 0; 0 255 0; 255 0 255; 200 170 50; 0 255 255; 120 50 200; 0 0 0]/255;       
        numbers=[2         151         300         449         598         747         896        1045        1194];
        leg='1234567';
        abs_min=0; abs_max=0;
        xlab=2.4:0.001:7.2;
        for j=1:9
            cnt=numbers(j);
            for i=0:20:130
                subplot(3,3,j)
                tmp=zeros(10,7740);
                kk=1;
                for k=cnt+i:2:cnt+i+19
                    tmp(kk,:)=convolved(spikes{k},40,7500);
                    kk=kk+1;
                end               
                tmp=mean(tmp,1);
                tmp=tmp(2400:7200)-mean(tmp(500:2000));
                plot(xlab,tmp,'color',colors(i/20+1,:))
                hold on
                if j>2
                    abs_min=min(abs_min,min(min(tmp)));
                    abs_max=max(abs_max,max(max(tmp)));
                end
            end
            title(['ND ',int2str(10-j)])
            xlabel('s')
            ylabel('a.u.')
            if j==1
                legend(leg')
            end
        end
        for j=3:9
            subplot(3,3,j)
            axis([2.4 7.2 abs_min abs_max])
        end
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(['20110818  CH',int2str(ch),', spikes minus min, boxes: ND 9-1, lines in boxes: mean of each 10 subs. trials of CONTRAST 2'],'FontSize',15,'FontWeight','bold')
        saveas(gcf,['F:\20110818\processing\pictures\spikes_by_channel_min\20110818_CH',int2str(ch),'.bmp']);
    end
end

% plot lfp by channel 27.05
for ch=1:60
    if ch~=15
        load(['F:\20110526\processing\preprocessed\CH',int2str(ch),'_0526'],'lfp');
        close all
        figure(1)
        set(gcf,'Position',[1  31 1680 946]);
        colors=[0 0 255; 255 0 0; 0 255 0; 255 0 255; 200 170 50; 0 255 255; 120 50 200; 0 0 0]/255;
        sf=250;%sampling frequency
        passL=30;orderL=10;
        fNormL = passL/(sf/2);
        [aL,bL] = butter(orderL, fNormL, 'low');
        numbers=[2         152         300         452         602         752         902        1052        1202];
        leg='1234567';
        abs_min=0; abs_max=0;
        xlab=2.4:0.004:7.2;
        for j=1:9
            cnt=numbers(j);
            for i=0:20:130
                subplot(3,3,j)
                tmp=filtfilt(aL, bL, mean(lfp(6000:10:18000,cnt+i:2:cnt+i+19),2));
                plot(xlab,tmp,'color',colors(i/20+1,:))
                hold on
                abs_min=min(abs_min,min(min(tmp)));
                abs_max=max(abs_max,max(max(tmp)));
            end
            title(['ND ',int2str(10-j)])
            xlabel('s')
            ylabel('a.u.')
            if j==1
                legend(leg')
            end
        end
        for j=1:9
            subplot(3,3,j)
            axis([2.4 7.2 abs_min abs_max])
        end
        % subplot('Position',[0.2 0.97 0.6 0.02])
        % set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        % text(0.02,0.5,['20110806  CH',int2str(ch),', lfp <30Hz, ND 9-1, '])
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(['20110526  CH',int2str(ch),', LFP, low pass filt. 30Hz, boxes: ND 9-1, lines in boxes: mean of each 10 subs. trials of CONTRAST 2'],'FontSize',15,'FontWeight','bold')
        saveas(gcf,['F:\20110526\processing\pictures\lfp_by_channel\20110526_CH',int2str(ch),'.bmp']);
    end
end

% plot spikes by channel minus mean spont 27.05
for ch=1:60
    if ch~=15
        load(['F:\20110526\processing\preprocessed\CH',int2str(ch),'_0526'],'spikes');
        close all
        figure(1)
        set(gcf,'Position',[1  31 1680 946]);
        colors=[0 0 255; 255 0 0; 0 255 0; 255 0 255; 200 170 50; 0 255 255; 120 50 200; 0 0 0]/255;       
        numbers=[2         152         300         452         602         752         902        1052        1202];
        leg='1234567';
        abs_min=0; abs_max=0;
        xlab=2.4:0.001:7.2;
        for j=1:9
            cnt=numbers(j);
            for i=0:20:130
                subplot(3,3,j)
                tmp=zeros(10,7740);
                kk=1;
                for k=cnt+i:2:cnt+i+19
                    tmp(kk,:)=convolved(spikes{k},40,7500);
                    kk=kk+1;
                end               
                tmp=mean(tmp,1);
                tmp=tmp(2400:7200)-mean(tmp(500:2000));
                plot(xlab,tmp,'color',colors(i/20+1,:))
                hold on
                if j>2
                    abs_min=min(abs_min,min(min(tmp)));
                    abs_max=max(abs_max,max(max(tmp)));
                end
            end
            title(['ND ',int2str(10-j)])
            xlabel('s')
            ylabel('a.u.')
            if j==1
                legend(leg')
            end
        end
        for j=3:9
            subplot(3,3,j)
            axis([2.4 7.2 abs_min abs_max])
        end
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(['20110526  CH',int2str(ch),', spikes minus min, boxes: ND 9-1, lines in boxes: mean of each 10 subs. trials of CONTRAST 2'],'FontSize',15,'FontWeight','bold')
        saveas(gcf,['F:\20110526\processing\pictures\spikes_by_channel_min\20110526_CH',int2str(ch),'.bmp']);
    end
end


%% Other stuff

% plots cleaned averaged lfp from the retina
% cleaning: remove "bad" channels; normalize amplitude by max achieved on the channel.
badCH=cell(1,8);
badCH{1}=[15 45];
badCH{2}=[15 11 45];
badCH{3}=[15 22 23 25 33 38];
badCH{4}=[15];
badCH{5}=[15 1 2 3 5 6 8 9 13 14 25 26 33 34 59 60];
badCH{6}=[15 30];
badCH{7}=[15 16 17 30 34 38 43 44 45 47 48 57 58];
badCH{8}=[15 3 4 5 6 7 8 9 27 43 47 48 53 57 58];
dates=['0524';'0525';'0526';'0527';'0805';'0806';'0809';'0818']
for cur=8
date=dates(cur,:);
ch_cnt=0;
tmp=zeros(1201,7,9,60-length(badCH{cur}));
sf=250;%sampling frequency
passL=30;orderL=10;
fNormL = passL/(sf/2);
[aL,bL] = butter(orderL, fNormL, 'low');
tic
for ch=1:60
    if isempty(find(badCH{cur}==ch, 1))
        ch_cnt=ch_cnt+1
        if cur<5
            load(['F:\2011',date,'\processing\preprocessed\CH',int2str(ch),'_',date],'lfp');
            numbers=[2         152         300         452         602         752         902        1052        1202];
        else
             load(['F:\2011',date,'\processing\preprocessed\lfp\quick_CH',int2str(ch),'_',date],'lfp');
             numbers=[2         151         300         449         598         747         896        1045        1194];
        end
        for j=1:9
            cnt=numbers(j);
            for i=0:20:130
                tmp(:,i/20+1,j,ch_cnt)=filtfilt(aL, bL, mean(lfp(6000:10:18000,cnt+i:2:cnt+i+19),2));
            end
        end
    end
end
toc
 % Now, normalization by 0-to-min distance on each channel separately
for i=1:size(tmp,4)
    tmptmp=tmp(:,:,:,i);
    a=min(min(min(tmptmp)));
    tmp(:,:,:,i)=-tmptmp/a;  
end
save(['F:\comparison_in_ND_lfp\',date],'tmp');
close(figure(1))
figure(1)
set(gcf,'Position',[1  31 1680 946]);
colors=[0 0 255; 255 0 0; 0 255 0; 255 0 255; 200 170 50; 0 255 255; 120 50 200; 0 0 0]/255;
leg='1234567';
xlab=2.4:0.004:7.2;
for j=1:9    
    for i=1:7
        subplot(3,3,j)
        a=reshape(tmp(:,i,j,:),1201,size(tmp,4));
        plot(xlab,mean(a,2),'color',colors(i,:))
        hold on
    end
    title(['ND ',int2str(10-j)])
    xlabel('s')
    ylabel('a.u.')
    if j==1
        legend(leg')
    end
end
for j=1:9
    subplot(3,3,j)
    axis tight
end
subplot('Position',[0.5 0.95 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(['2011',date,'  clean channels averaged, LFP, low pass filt. 30Hz, boxes: ND 9-1, lines in boxes: mean of each 10 subs. trials of CONTRAST 2'],'FontSize',15,'FontWeight','bold')
saveas(gcf,['F:\comparison_in_ND_lfp\tight_axis\tight_2011',date,'.bmp']);
close(figure(2))
figure(2)
set(gcf,'Position',[1  31 1680 946]);
colors=[0 0 255; 255 0 0; 0 255 0; 255 0 255; 200 170 50; 0 255 255; 120 50 200; 0 0 0]/255;
leg='1234567';
xlab=2.4:0.004:7.2;
for j=1:9    
    for i=1:7
        subplot(3,3,j)
        a=reshape(tmp(:,i,j,:),1201,size(tmp,4));
        plot(xlab,mean(a,2),'color',colors(i,:))
        hold on
    end
    title(['ND ',int2str(10-j)])
    xlabel('s')
    ylabel('a.u.')
    if j==1
        legend(leg')
    end
end
for j=1:9
    subplot(3,3,j)
    axis([2.4 7.2 -1 max(max(max(max(tmp))))])
end
subplot('Position',[0.5 0.95 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(['2011',date,'  clean channels averaged, LFP, low pass filt. 30Hz, boxes: ND 9-1, lines in boxes: mean of each 10 subs. trials of CONTRAST 2'],'FontSize',15,'FontWeight','bold')
saveas(gcf,['F:\comparison_in_ND_lfp\equal_axis\equal_2011',date,'.bmp']);

end


%%%% END %%%%%



% plots cleaned averaged SPIKES (convolved) from the retina
% cleaning: remove "bad" channels. tmp(:,i/20+1,j,ch_cnt): time count,
% trial set, ND, channel (only clean ones)
badCH=cell(1,8);
badCH{1}=[15 5 7 9 10 12 35 36 37 39 40 42];
badCH{2}=[15 33];
badCH{3}=[15 24 25 38];
badCH{4}=[15 12];
badCH{5}=[15 1 2 6 8 26 59 60];
badCH{6}=[15];
badCH{7}=[15 43];
badCH{8}=[15 1 3 4 7 24 26 27 28 30 33 34 38 43 44 45 46 47 48 53 57 58];
dates=['0524';'0525';'0526';'0527';'0805';'0806';'0809';'0818']
for cur=8
date=dates(cur,:);
ch_cnt=0;
tmp=zeros(4801,7,9,60-length(badCH{cur}));
tic
for ch=1:60
    if isempty(find(badCH{cur}==ch, 1))
        ch_cnt=ch_cnt+1
        if cur<5
            load(['F:\2011',date,'\processing\preprocessed\CH',int2str(ch),'_',date],'spikes');
            numbers=[2         152         300         452         602         752         902        1052        1202];
        else
             load(['F:\2011',date,'\processing\preprocessed\spikes\quick_CH',int2str(ch),'_',date],'spikes');
             numbers=[2         151         300         449         598         747         896        1045        1194];
        end
        for j=1:9
            cnt=numbers(j);
            for i=0:20:130
                kk=1;
                tmp1=zeros(10,7740);
                for k=cnt+i:2:cnt+i+19
                    tmp1(kk,:)=convolved(spikes{k},40,7500);
                    kk=kk+1;
                end               
                tmp1=mean(tmp1,1);
                tmp(:,i/20+1,j,ch_cnt)=tmp1(2400:7200)-mean(tmp1(500:2000));
            end
        end
    end
end
toc

save(['F:\comparison_in_ND_spikes\',date],'tmp');
close(figure(1))
figure(1)
set(gcf,'Position',[1  31 1680 946]);
colors=[0 0 255; 255 0 0; 0 255 0; 255 0 255; 200 170 50; 0 255 255; 120 50 200; 0 0 0]/255;
leg='1234567';
xlab=2.4:0.001:7.2;
for j=1:9    
    for i=1:7
        subplot(3,3,j)
        a=reshape(tmp(:,i,j,:),4801,size(tmp,4));
        plot(xlab,mean(a,2),'color',colors(i,:))
        hold on
    end
    title(['ND ',int2str(10-j)])
    xlabel('s')
    ylabel('a.u.')
    if j==1
        legend(leg')
    end
end
for j=1:9
    subplot(3,3,j)
    axis tight
end
subplot('Position',[0.5 0.95 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(['2011',date,'  clean channels averaged, spikes convolved; boxes: ND 9-1, lines in boxes: mean of each 10 subs. trials of CONTRAST 2'],'FontSize',15,'FontWeight','bold')
saveas(gcf,['F:\comparison_in_ND_spikes\tight\tight_2011',date,'.bmp']);
close(figure(2))
figure(2)
set(gcf,'Position',[1  31 1680 946]);
colors=[0 0 255; 255 0 0; 0 255 0; 255 0 255; 200 170 50; 0 255 255; 120 50 200; 0 0 0]/255;
leg='1234567';
xlab=2.4:0.001:7.2;
m=0;
for j=1:9    
    for i=1:7
        subplot(3,3,j)
        a=reshape(tmp(:,i,j,:),4801,size(tmp,4));
        plot(xlab,mean(a,2),'color',colors(i,:))
        m=max(max(mean(a,2)),m);
        hold on
    end
    title(['ND ',int2str(10-j)])
    xlabel('s')
    ylabel('a.u.')
    if j==1
        legend(leg')
    end
end
for j=1:9
    subplot(3,3,j)
    axis([2.4 7.2 0 m])
end
subplot('Position',[0.5 0.95 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(['2011',date,'  clean channels averaged, spikes convolved; boxes: ND 9-1, lines in boxes: mean of each 10 subs. trials of CONTRAST 2'],'FontSize',15,'FontWeight','bold')
saveas(gcf,['F:\comparison_in_ND_spikes\equal\equal_2011',date,'.bmp']);

end



%%%%% END %%%%%



dates=['0524';'0525';'0526';'0527';'0805';'0806';'0809']
for cur=3:7
date=dates(cur,:);
load(['F:\comparison_in_ND_spikes\',date],'tmp');

close(figure(2))
figure(2)
set(gcf,'Position',[1  31 1680 946]);
m=1;n=10000;
for j=1:9
    for i=1:7
        a=mean(tmp(:,i,j,:),4);
        onmax(i)=max(a(500:1500));
        offmax(i)=max(a(2500:3500));
    end
    subplot(3,3,j)
    plot(onmax)
    hold on
    plot(offmax,'color','r')
    title(['ND ',int2str(10-j)])
    m=max(max(onmax(i)),m);
    m=max(max(offmax(i)),m);
    n=min(min(onmax(i)),n);
    n=min(min(offmax(i)),n);
end
for j=1:9
    subplot(3,3,j)
    axis([1 7 n*1.1 m*1.1])
end
subplot('Position',[0.5 0.95 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(['2011',date,'  on (blue) and off(red) spiking reaction maximas; boxes: ND 9-1, lines in boxes: mean of each 10 subs. trials of CONTRAST 2'],'FontSize',15,'FontWeight','bold')
saveas(gcf,['F:\comparison_in_ND_spikes\maxima\',date,'.bmp'])
end



close(figure(1))
figure(1)
set(gcf,'Position',[1  31 1680 946]);

dates=['0524';'0525';'0526';'0527';'0805';'0806';'0809']
for cur=3:7
date=dates(cur,:);
load(['F:\comparison_in_ND_spikes\',date],'tmp');


m=1;n=10000;
for j=1:9
    for i=1:7
        a=mean(tmp(:,i,j,:),4);
        onmax(i)=max(a(500:1500));
        offmax(i)=max(a(2500:3500));
    end
    subplot(3,3,j)
%     plot(onmax-min(onmax))
    hold on
    plot(offmax-min(offmax),'color','r')
    title(['ND ',int2str(10-j)])
%     m=max(max(onmax(i)),m);
%     m=max(max(offmax(i)),m);
%     n=min(min(onmax(i)),n);
%     n=min(min(offmax(i)),n);
end
% for j=9
%     subplot(3,3,j)
%     axis([1 7 n*1.1 m*1.1])
% end
end
subplot('Position',[0.5 0.95 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(['2011',date,'  on (blue) and off(red) spiking reaction maximas; boxes: ND 9-1, lines in boxes: mean of each 10 subs. trials of CONTRAST 2'],'FontSize',15,'FontWeight','bold')
saveas(gcf,['F:\comparison_in_ND_spikes\maxima\',date,'.bmp'])









%plots raw data from date 'date', channel 'ch', listed in 'list', and
%overlay it with convolved spiking rate
figure(4)
ch='22'; date='0805';
list=[350 460 580 850] % ND7, ND6 begin, ND6 end, ND3
stim=[2984.6, 4968.7];
load(['F:\2011',date,'\processing\preprocessed\spikes\quick_CH',ch,'_',date])
for i=1:length(list)
    b=list(i);
    if b<999
        k=['0',int2str(b)];
    else
        k=int2str(b);
    end
    subplot(length(list),1,i);

    c=convolved(spikes{b},40,0);
    c=c-mean(c(1000:2000));
    c=c(120:end-120);
    c(1)=-1000;
    c(end)=1000;
    
    load(['F:\2011',date,'\MEA_bin\channels\channel',ch,'\2011',date,'_CH',ch,'_',k,'_inND_.mat'])
    a=filtfilt(aH,bH,data);
    a=a(1:length(c)*25);
    xlab=1:length(a);
    xlab=xlab/(length(a)/length(c));
    
    [AX,H1,H2]=plotyy(xlab,a,1:length(c),c);
    set(H1,'color','b');
    set(H2,'color','r','LineWidth',2.5);
%     axis([0 xlab(end) -750 750 ])
%     line([0,xlab(end)],[-500,-500],'color','k')
%     line([0,xlab(end)],[500,500],'color','k')
%     line([stim(1),stim(1)],[-500, 500],'color','m','LineWidth',3)
%     line([stim(2),stim(2)],[-500, 500],'color','m','LineWidth',3)
    title(k)    
    
end






figure(4)
load('F:\20110806\MEA_bin\channels\channel46\20110806_CH46_0286_inND_.mat')
a=filtfilt(aH,bH,data);
dataHPF=a.*a.*a;
dataHPF=dataHPF';
spikesTotal=maxSpikesAllowed+1;
coef=1;
while spikesTotal>maxSpikesAllowed
    test=dataHPF;
    test1=test;
    test1(test>-coef*std(test([1:10000 end-10000:end])))=0;
    test=diff(test1);
    fr=[0; (test==test1(2:end) & test)];
    spikesTotal=nnz(fr);
    coef=coef+0.5;
end
fr=fr(1:len*25);
% get spike times
fr=find(fr);
subplot(2,1,1)
tmp=0;
hold on
for i=1:length(fr)
%     tmp=min(a(fr(i)-20:fr(i)+20));
    plot(a(fr(i)-20:fr(i)+20)-a(fr(i)));
end
subplot(2,1,2)
hold on
for i=1:length(fr)
%     tmp=min(data(fr(i)-20:fr(i)+20));
    plot(data(fr(i)-20:fr(i)+20)-data(fr(i)));
end


