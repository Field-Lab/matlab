%% PLotting peaks averaged from all channels
% change "parameterToPlot" to plot: 1) amplitude of peaks; 2)mean of std of peaks
% 3) latencies of peaks


%% all channels mean


stim=[2951.2, 4935.3];
stimLFP=floor(stim);
stimBin=floor(stim);
m='9876543211234543210'; % ND list
start2plot=1;
filesToPlot=length(m)*148;

contrast=0;

if contrast==1
    contr=2:2:filesToPlot;
    sk=74;
    cont='ON';
elseif contrast==2
    contr=1:2:filesToPlot;
    sk=74;    
    cont='OFF';
else       
    contr=1:filesToPlot;
    sk=148;
    cont='ON and OFF';
end

parameterToPlot=1; % 1 - amplitudes; 2 - std; 3 - LP

OnSpikesAll=[];
OffSpikesAll=[];
OnLFPAll=[];
OffLFPminimaAll=[];
OffLFPmaximaAll=[];
spontSpikesAll=[];
ch_list=[2 5:10 12 13 16 19:26 28:37 39:42 44 46 49 51 52 54 56 59];
for i=1:60
    if ~isempty(find(ch_list==i, 1))
        load(['F:\20110913\processing\ready_all\ready_CH',int2str(i),'_0913.mat']);
        OnSpikesAll=[OnSpikesAll; OnSpikes(parameterToPlot,start2plot:filesToPlot)];
        OffSpikesAll=[OffSpikesAll; OffSpikes(parameterToPlot,start2plot:filesToPlot)];
        OnLFPAll=[OnLFPAll; OnLFP(parameterToPlot,start2plot:filesToPlot)];
        OffLFPminimaAll=[OffLFPminimaAll; OffLFPminima(parameterToPlot,start2plot:filesToPlot)];
        OffLFPmaximaAll=[OffLFPmaximaAll; OffLFPmaxima(parameterToPlot,start2plot:filesToPlot)];
        if parameterToPlot~=3
            spontSpikesAll=[spontSpikesAll; spontSpikes(parameterToPlot,start2plot:filesToPlot)];
        else
            spontSpikesAll=[spontSpikesAll; zeros(1,filesToPlot)];
        end
    end
end

for i=1:length(ch_list)
    OnSpikesAll(i,:)=OnSpikesAll(i,:)-spontSpikesAll(i,:);
    OffSpikesAll(i,:)=OffSpikesAll(i,:)-spontSpikesAll(i,:);
    OnSpikesAll(i,:)=(OnSpikesAll(i,:)-min(OnSpikesAll(i,:)))/max(OnSpikesAll(i,:)-min(OnSpikesAll(i,:)));
    OffSpikesAll(i,:)=(OffSpikesAll(i,:)-min(OffSpikesAll(i,:)))/max(OffSpikesAll(i,:)-min(OffSpikesAll(i,:)));
    OnLFPAll(i,:)=(OnLFPAll(i,:)-min(OnLFPAll(i,:)))/max(OnLFPAll(i,:)-min(OnLFPAll(i,:)));
    OffLFPminimaAll(i,:)=(OffLFPminimaAll(i,:)-min(OffLFPminimaAll(i,:)))/max(OffLFPminimaAll(i,:)-min(OffLFPminimaAll(i,:)));
    OffLFPmaximaAll(i,:)=(OffLFPmaximaAll(i,:)-min(OffLFPmaximaAll(i,:)))/max(OffLFPmaximaAll(i,:)-min(OffLFPmaximaAll(i,:)));
    spontSpikesAll(i,:)=(spontSpikesAll(i,:)-min(spontSpikesAll(i,:)))/max(spontSpikesAll(i,:)-min(spontSpikesAll(i,:)));
end

titles=cell(1,6);
titles{1}='Maxima of Spiking On-reaction, spont substracted';
titles{2}='Maxima of Spiking Off-reaction, spont substracted';
titles{3}='Minima of LFP On-reaction';
titles{4}='Minima of LFP Off-reaction';
titles{5}='Maxima of LFP Off-reaction';
titles{6}='Spontaneous spikes';

CoF=[0.5,0.5,0.5];
OnSpikes=mean(OnSpikesAll,1);
OffSpikes=mean(OffSpikesAll,1);
OnLFP=mean(OnLFPAll,1);
OffLFPminima=mean(OffLFPminimaAll,1);
OffLFPmaxima=mean(OffLFPmaximaAll,1);
spontSpikes=mean(spontSpikesAll,1);

switch parameterToPlot
    case 1
        toPlot=[OnSpikes(contr);OffSpikes(contr);1-OnLFP(contr);1-OffLFPminima(contr);OffLFPmaxima(contr);spontSpikes(contr)];
    case 2
        toPlot=[OnSpikes;OffSpikes;OnLFP;OffLFPminima;OffLFPmaxima;spontSpikes];
    case 3
        toPlot=[OnSpikes*20;OffSpikes*20;OnLFP*0.4;OffLFPminima*0.4;OffLFPmaxima*0.4];
end

a=size(toPlot,1);
figure(1)
set(gcf,'position',[10 50 1660 900])
for j=1:a
    subplot(a,1,j)
    plot(toPlot(j,:))
    fill(0:size(toPlot,2)+1,[0,toPlot(j,:),0],CoF);
    for i=1:length(m)
        text((i-1)*sk+74/sk,1.1,m(i),'color','r','FontSize',11);
        line([(i-1)*sk+1,(i-1)*sk+1],[0,1],'color','c','LineWidth',2);
    end
    axis tight;
    title(['20110913 Contrast ',cont,' ', titles{j}])
end

saveas(gcf,['F:\20110913\processing\pictures\all\all_',cont,'_0913.bmp']);
close (gcf)



%% Plot by Channel
cd('F:\20110913\processing')
ch_list=[2 5:10 12 13 16 19:26 28:37 39:42 44 46 49 51 52 54 56 59];
% plot lfp by channel
for ch=1:60
    if ~isempty(find(ch_list==ch, 1))
        load(['F:\20110913\processing\preprocessed\lfp\quick_1_1350_CH',int2str(ch),'_0913']);
        data=lfp;
        load(['F:\20110913\processing\preprocessed\lfp\back_1355_2105_CH',int2str(ch),'_0913']);
        tmp=lfp;
        lfp(:,1:2:end)=lfp(:,2:2:end);
        lfp(:,2:2:end)=tmp(:,1:2:end); % change so that the contrasts are the same order
        data=[data lfp];
        load(['F:\20110913\processing\preprocessed\lfp\quick_2255_3004_CH',int2str(ch),'_0913']);
        data=[data lfp];
        clear lfp tmp
        close all
        figure(1)
        set(gcf,'Position',[1  31 1680 946]);
        colors=[0 0 255; 255 0 0; 0 255 0; 255 0 255; 200 170 50; 0 255 255; 120 50 200; 0 0 0]/255;
        sf=1000;%sampling frequency
        passL=30;orderL=10;
        fNormL = passL/(sf/2);
        [aL,bL] = butter(orderL, fNormL, 'low');
        numbers=[148:148:1184 1776 1628 1480 1332 2072:148:2516];        
        leg='1234567';
        tit='8765432143214321';
        abs_min=0; abs_max=0;
        xlab=2.5:0.001:7;
        for j=1:16
            cnt=numbers(j);
            for i=0:20:130
                subplot(4,4,j)
                tmp=filtfilt(aL, bL, mean(data(:,cnt+i:2:cnt+i+19),2));
                tmp=tmp(2500:7000);
                plot(xlab,tmp,'color',colors(i/20+1,:))
                hold on
                abs_min=min(abs_min,min(min(tmp)));
                abs_max=max(abs_max,max(max(tmp)));
            end
            title(['ND ',tit(j)])
            xlabel('s')
            ylabel('a.u.')
            if j==1
                legend(leg')
            end
        end
        for j=1:16
            subplot(4,4,j)
            axis([2.5 7 abs_min abs_max])
        end
        % subplot('Position',[0.2 0.97 0.6 0.02])
        % set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        % text(0.02,0.5,['20110913  CH',int2str(ch),', lfp <30Hz, ND 9-1, '])
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(['20110913  CH',int2str(ch),', LFP, low pass filt. 30Hz, boxes: ND 8-1-1-4-4-1, lines in boxes: mean of each 10 subs. trials, BG25, ST20'],'FontSize',15,'FontWeight','bold')
        saveas(gcf,['F:\20110913\processing\pictures\lfp_by_channel_all\BG25ST20_20110913_CH',int2str(ch),'.bmp']);
    end
end
% plot spikes by channel minus mean spont
for ch=1:60
    if ch~=15
        load(['F:\20110913\processing\preprocessed\spikes\back_1355_2105_CH',int2str(ch),'_0913']);
        close all
        figure(1)
        set(gcf,'Position',[1  31 1680 946]);
        colors=[0 0 255; 255 0 0; 0 255 0; 255 0 255; 200 170 50; 0 255 255; 120 50 200; 0 0 0]/255;       
        numbers=2:148:750;
        leg='1234567';
        abs_min=0; abs_max=0;
        xlab=2.4:0.001:7.2;
        for j=1:5
            cnt=numbers(j);
            for i=0:20:130
                subplot(2,3,j)
                tmp=zeros(10,7940);
                kk=1;
                for k=cnt+i:2:cnt+i+19
                    tmp(kk,:)=convolved(spikes{k},40,7700);
                    kk=kk+1;
                end               
                tmp=mean(tmp,1);
                tmp=tmp(2400:7200)-mean(tmp(500:2000));
                plot(xlab,tmp,'color',colors(i/20+1,:))
                hold on
                if j>0
                    abs_min=min(abs_min,min(min(tmp)));
                    abs_max=max(abs_max,max(max(tmp)));
                end
            end
            title(['ND ',int2str(j)])
            xlabel('s')
            ylabel('a.u.')
            if j==1
                legend(leg')
            end
        end
        for j=1:5
            subplot(2,3,j)
            axis([2.4 7.2 abs_min abs_max])
        end
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(['20110913  CH',int2str(ch),', spikes minus min, boxes: ND 1-5, lines in boxes: mean of each 10 subs. trials of CONTRAST 1'],'FontSize',15,'FontWeight','bold')
        saveas(gcf,['F:\20110913\processing\pictures\spikes_by_channel_min_back\Contrast1_20110913_CH',int2str(ch),'.bmp']);
    end
end






