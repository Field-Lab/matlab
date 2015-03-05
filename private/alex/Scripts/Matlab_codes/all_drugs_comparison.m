%% for comparison analysis
% Mean and std of ampl and LP in control, drug, wash - always first ND4 (ND3 for 0319).
% Significance - Friedman's test.
% procedure. 1) separate for each channel: normalize amplitude and latency by corresponding mean of control,
% which is taken for 100%; calculate mean and std for all conditions, each
% value shows therefore the mean and variability of peaks on each channel.
% Calculate Friedman test for each channel - matrix for all conditions.
% 2) calculate mean of all values all channels. "Std" reflects then the average
% variablity of amplitude/latency across all 59 channels for each
% condition, "Fried" reflects the average significance across all 59
% channels.
% structure of "allValues" (in order needed for drawning)
% 1 row: mean of OnLFP amplitude
% 2 row: std of OnLFP amplitude
% 3 row: mean of OffLPFminima amplitude
% 4 row: std of OffLPFminima amplitude
% 5 row: mean of OffLPFmaxima amplitude
% 6 row: std of OffLPFmaxima amplitude
% 7 row: mean of OnLFP latency
% 8 row: std of OnLFP latency
% 9 row: mean of OffLPFminima latency
% 10 row: std of OffLPFminima latency
% 11 row: mean of OffLPFmaxima latency
% 12 row: std of OffLPFmaxima latency
% structure of variable "fried"
% rows:
% 1 row:  OnLFP amplitude
% 2 row:  OnLFP latency
% 3 row:  OffLPFminima amplitude
% 4 row:  OffLPFminima latency
% 5 row:  OffLPFmaxima amplitude
% 6 row:  OffLPFmaxima latency
% 7 row:  index of condition 1 in pair compared in this column
% 8 row:  index of condition 2 in pair compared in this column
% columns: vectorized matrix of pairwise conditions comparison along 1, then 2 etc. -
% [(1,2),(1,3),...(1,end),(2,3),..(end-1),(end)].

dates=['0308';'0311';'0315';'0316';'0319'];

for dateCNT=1:5
    dateR=dates(dateCNT,:);

    switch dateR
        case '0319'
            timesStart=[50 87 127 227];
            timesEnd=[63 100 142 242];
            titles=cell(1,4);
            titles{1}='Control';
            titles{2}='Cocktail 20';
            titles{3}='Cocktail 100';
            titles{4}='Wash';
            fileNumbers='561-849';
            artifacts=[3,19,20,24,55,65,97,98,103,126,129,137,138,143,161,246,247,248,249,265,266];
        case '0308'
            timesStart=[3 30 93];
            timesEnd=[14 45 103];
            titles=cell(1,3);
            titles{1}='Control';
            titles{2}='APB 10';
            titles{3}='Wash';
            fileNumbers='201-633';
            artifacts=[208,209,215,216,290,291,292,293,300,302,303,304]-200;
        case '0311'
            timesStart=[3 29 137];
            timesEnd=[11 43 159];
            titles=cell(1,3);
            titles{1}='Control';
            titles{2}='APB 10';
            titles{3}='Wash';
            fileNumbers='201-462';
            artifacts=[207,208,209,215,216,224,244,260,275,286,300,380,441,455]-199;
        case '0315'
            timesStart=[3 50 234];
            timesEnd=[18 67 260];
            titles=cell(1,3);
            titles{1}='Control';
            titles{2}='CPP/NBQX 10';
            titles{3}='Wash';
            fileNumbers='236-723';
            artifacts=[207,216,218,224,246,257,263,273,276,305,326,369,464,521]-200;
        case '0316'
            timesStart=[100 150 300 460 502 600];
            timesEnd=[123 168 313 475 513 612];
            titles=cell(1,6);
            titles{1}='Control';
            titles{2}='CPP/NBQX 10';
            titles{3}='Wash';
            titles{4}='CPP/NBQX 20';
            titles{5}='Cocktail 20';
            titles{6}='Wash';
            fileNumbers='433-1132';
            artifacts=[121,122,265,266,365,393,394,443,533];
    end
    parName=cell(1,3);
    parName{1}='On LFP';
    parName{2}='Off LFP negative';
    parName{3}='Off LFP positive';
    
    
    allValues=zeros(12,length(timesStart),59);
    fried=zeros(6,length(timesStart),59);
    cnt=1;
    for l=1:60
        if l~=15
            Channel=int2str(l);
            %         display(Channel)
            fileName=['F:\2011',dateR,'\processing\ready\',fileNumbers,'_CH',Channel,'_',dateR,'.ready.mat'];
            load(fileName);
            %artifacts exclusion
            for j=1:length(artifacts)
                k=artifacts(j);
                OnLFP(:,k)=OnLFP(:,k-1);
                OffLFPminima(:,k)=OffLFPminima(:,k-1);
                OffLFPmaxima(:,k)=OffLFPmaxima(:,k-1);
            end
            % normalization to 100 % of control
            normOnLFP=mean(OnLFP(1,timesStart(1):timesEnd(1)));
            OnLFP(1,:)=OnLFP(1,:)/normOnLFP*100;
            normOffLFPminima=mean(OffLFPminima(1,timesStart(1):timesEnd(1)));
            OffLFPminima(1,:)=OffLFPminima(1,:)/normOffLFPminima*100;
            normOffLFPmaxima=mean(OffLFPmaxima(1,timesStart(1):timesEnd(1)));
            OffLFPmaxima(1,:)=OffLFPmaxima(1,:)/normOffLFPmaxima*100;
            
            normOnLFP=mean(OnLFP(3,timesStart(1):timesEnd(1)));
            OnLFP(3,:)=OnLFP(3,:)/normOnLFP*100;
            normOffLFPminima=mean(OffLFPminima(3,timesStart(1):timesEnd(1)));
            OffLFPminima(3,:)=OffLFPminima(3,:)/normOffLFPminima*100;
            normOffLFPmaxima=mean(OffLFPmaxima(3,timesStart(1):timesEnd(1)));
            OffLFPmaxima(3,:)=OffLFPmaxima(3,:)/normOffLFPmaxima*100;
            
            for i=1:length(timesStart)
                allValues(1,i,cnt)=mean(OnLFP(1,timesStart(i):timesEnd(i)));
                allValues(2,i,cnt)=std(OnLFP(1,timesStart(i):timesEnd(i)));
                allValues(3,i,cnt)=mean(OffLFPminima(1,timesStart(i):timesEnd(i)));
                allValues(4,i,cnt)=std(OffLFPminima(1,timesStart(i):timesEnd(i)));
                allValues(5,i,cnt)=mean(OffLFPmaxima(1,timesStart(i):timesEnd(i)));
                allValues(6,i,cnt)=std(OffLFPmaxima(1,timesStart(i):timesEnd(i)));
                allValues(7,i,cnt)=mean(OnLFP(3,timesStart(i):timesEnd(i)));
                allValues(8,i,cnt)=std(OnLFP(3,timesStart(i):timesEnd(i)));
                allValues(9,i,cnt)=mean(OffLFPminima(3,timesStart(i):timesEnd(i)));
                allValues(10,i,cnt)=std(OffLFPminima(3,timesStart(i):timesEnd(i)));
                allValues(11,i,cnt)=mean(OffLFPmaxima(3,timesStart(i):timesEnd(i)));
                allValues(12,i,cnt)=std(OffLFPmaxima(3,timesStart(i):timesEnd(i)));
                
            end
            cnt1=1;
            for i=1:length(timesStart)-1
                for j=i+1:length(timesStart)
                    a=min(length(timesStart(i):timesEnd(i)),length(timesStart(j):timesEnd(j)));
                    fried(1,cnt1,cnt)=friedman([OnLFP(1,timesStart(i):timesStart(i)+a-1)',OnLFP(1,timesStart(j):timesStart(j)+a-1)'],1,'off');
                    fried(2,cnt1,cnt)=friedman([OnLFP(3,timesStart(i):timesStart(i)+a-1)',OnLFP(3,timesStart(j):timesStart(j)+a-1)'],1,'off');
                    fried(3,cnt1,cnt)=friedman([OffLFPminima(1,timesStart(i):timesStart(i)+a-1)',OffLFPminima(1,timesStart(j):timesStart(j)+a-1)'],1,'off');
                    fried(4,cnt1,cnt)=friedman([OffLFPminima(3,timesStart(i):timesStart(i)+a-1)',OffLFPminima(3,timesStart(j):timesStart(j)+a-1)'],1,'off');
                    fried(5,cnt1,cnt)=friedman([OffLFPmaxima(1,timesStart(i):timesStart(i)+a-1)',OffLFPmaxima(1,timesStart(j):timesStart(j)+a-1)'],1,'off');
                    fried(6,cnt1,cnt)=friedman([OffLFPmaxima(3,timesStart(i):timesStart(i)+a-1)',OffLFPmaxima(3,timesStart(j):timesStart(j)+a-1)'],1,'off');
                    fried(7,cnt1,cnt)=i;
                    fried(8,cnt1,cnt)=j;
                    cnt1=cnt1+1;
                end
            end
            cnt=cnt+1;
        end
    end
    allValues=mean(allValues,3);
    fried=mean(fried,3);
    
    for figNumber=1:3
        figure(figNumber)
        set(gcf,'Position',[368 272 1028 679])
        subplot('Position',[0.1 0.6 0.6 0.35])
        bar(allValues(1+(figNumber-1)*4,:),'FaceColor',[0.9 0.9 0.9])
        hold on
        errorbar(allValues(1+(figNumber-1)*4,:),allValues(2+(figNumber-1)*4,:),'LineStyle','none','color','k','LineWidth',2)
        hold off
        set(gca,'XTickLabel',titles)
        title(['2011-',dateR(1:2),'-',dateR(3:4),' Amplitude of ',parName{figNumber},', mean of control - 100%'])
        subplot('Position',[0.75 0.6 0.22 0.35])
        set(gca,'XTickLabel','','XTick',0,'YTickLabel','','YTick',0)
        [rows,colms]=find(fried(1:6,:)<0.01);
        tmp=colms(rows==(2*figNumber-1));
        sigLevel2plot='0.01';
        if isempty(tmp)
            [rows,colms]=find(fried(1:6,:)<0.05);
            tmp=colms(rows==(2*figNumber-1));
            sigLevel2plot='0.05';
        end            
        text(0.1,0.9,['Significance < ',sigLevel2plot],'fontweight','b','color','b')
        for i=1:length(tmp)
            a=fried(7,tmp(i));
            b=fried(8,tmp(i));
            text(0.1,0.8-0.1*i,[titles{a},'-',titles{b}])
        end
        
        subplot('Position',[0.1 0.1 0.6 0.35])
        bar(allValues(3+(figNumber-1)*4,:),'FaceColor',[0.9 0.9 0.9])
        hold on
        errorbar(allValues(3+(figNumber-1)*4,:),allValues(4+(figNumber-1)*4,:),'LineStyle','none','color','k','LineWidth',2)
        hold off
        set(gca,'XTickLabel',titles)
        title(['2011-',dateR(1:2),'-',dateR(3:4),' Latency of ',parName{figNumber},', mean of control - 100%'])
        subplot('Position',[0.75 0.1 0.22 0.35])
        set(gca,'XTickLabel','','XTick',0,'YTickLabel','','YTick',0)
        [rows,colms]=find(fried(1:6,:)<0.01);
        tmp=colms(rows==(2*figNumber));
        sigLevel2plot='0.01';
        if isempty(tmp)
            [rows,colms]=find(fried(1:6,:)<0.05);
            tmp=colms(rows==(2*figNumber));
            sigLevel2plot='0.05';
        end            
        text(0.1,0.9,['Significance < ',sigLevel2plot],'fontweight','b','color','b')
        for i=1:length(tmp)
            a=fried(7,tmp(i));
            b=fried(8,tmp(i));
            text(0.1,0.8-0.1*i,[titles{a},'-',titles{b}])
        end
        
        saveas(gcf,['F:\all_drugs_pictures\',dateR,'_',parName{figNumber},'.bmp']);
    end
    close all
    
end

%%%% code for if showSig %%%%
% first plot: 
% if showSig
%     tmp=colms(rows==(2*figNumber-1));
%     for i=1:length(tmp)
%         a=fried(7,tmp(i));
%         b=fried(8,tmp(i));
%         c=min([allValues(1+(figNumber-1)*4,a)-allValues(2+(figNumber-1)*4,a);allValues(1+(figNumber-1)*4,b)-allValues(2+(figNumber-1)*4,b)]);
%         line([a+0.1,b-0.1],[0.8*c,0.8*c],'color','r')
%     end
% end
% second plot:
% if showSig
%     tmp=colms(rows==(2*figNumber));
%     for i=1:length(tmp)
%         a=fried(7,tmp(i));
%         b=fried(8,tmp(i));
%         c=min([allValues(3+(figNumber-1)*4,a)-allValues(4+(figNumber-1)*4,a);allValues(3+(figNumber-1)*4,b)-allValues(4+(figNumber-1)*4,b)]);
%         line([a+0.1,b-0.1],[0.8*c,0.8*c],'color','r')
%     end
% end