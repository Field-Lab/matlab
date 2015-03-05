clear
date='20120902_1'
codeWord='ffflicker'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list')

FiringRate=zeros(310000,size(SpikeCount,1),size(SpikeCount,2));

for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,310000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
        end
    end
end

onOff=[1 -1 -1 1 -1 1 -1 -1 1 0 1 -1 -1 1 -1 1 1 0 1 -1 1 -1 0 0 -1 0 1 -1 1];

path2savePics=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/FiringRateAndLinearFilter/'];
if ~exist(path2savePics,'dir')
    mkdir(path2savePics);
end
figure
set(gcf,'position',[560          17        1677         931])
for cnt=1:size(FiringRate,3)
    if onOff(cnt)~=0
        for i=1:12
            subplot(3,4,i)
            hold off
            axis([0 60000 -80 120])
            if onOff(cnt)>0
                col='r';
            else
                col='b';
            end
            a=zeros(60000,5);
            p=1;
            for j=1:3600:18000
                tmp=FiringRate(correctedProtocols(j,1,1):correctedProtocols(j+3600,1,1),i,cnt);
                a(:,p)=tmp(1:60000);
                p=p+1;
            end
            a=mean(a,2);
            plot(a,col)
            hold on
            
            t=mean(LinearFilter(:,1:2:end,i,cnt),2);
            t=t-mean(t);
            t=t./sum(abs(t));
            ko=std(a(1:30000))/std(t);
            t=t*ko;
            tm=max(t);
            t=t-tm;
            plot(1:60:30000,t,col)
            
            t=mean(LinearFilter(:,2:2:end,i,cnt),2);
            t=t-mean(t);
            t=t./sum(abs(t));
            t=t*ko;
            t=t-tm;
            plot(30001:60:60000,t,'k')
        end
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(names{cnt},'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2savePics,names{cnt},'.png']);
    end
end



%% 20120902_2
clear
date='20120902_2'
codeWord='ffflicker'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

FiringRate=zeros(310000,size(SpikeCount,1),size(SpikeCount,2));

for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,310000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
        end
    end
end
% for i=1:24
%     subplot(3,8,i)
%     plot(mean(LinearFilter(:,1:2:end,3,i),2))
% end
onOff=[1 1 -1 0 0 0 1 1 1 -1 -1 -1 0 0 1 1 -1 1 0 -1 1 -1 1 1];

path2savePics=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/FiringRateAndLinearFilter/'];
if ~exist(path2savePics,'dir')
    mkdir(path2savePics);
end
figure
set(gcf,'position',[560          17        1677         931])
for cnt=1:size(FiringRate,3)
    if onOff(cnt)~=0
        for i=1:12
            subplot(3,4,i)
            hold off
            axis([0 60000 -80 120])
            if onOff(cnt)>0
                col='r';
            else
                col='b';
            end
            a=zeros(60000,5);
            p=1;
            for j=1:3600:18000
                tmp=FiringRate(correctedProtocols(j,1,1):correctedProtocols(j+3600,1,1),i,cnt);
                a(:,p)=tmp(1:60000);
                p=p+1;
            end
            a=mean(a,2);
            plot(a,col)
            hold on
            
            t=mean(LinearFilter(:,1:2:end,i,cnt),2);
            t=t-mean(t);
            t=t./sum(abs(t));
            ko=std(a(1:30000))/std(t);
            t=t*ko;
            tm=max(t);
            t=t-tm;
            plot(1:60:30000,t,col)
            
            t=mean(LinearFilter(:,2:2:end,i,cnt),2);
            t=t-mean(t);
            t=t./sum(abs(t));
            t=t*ko;
            t=t-tm;
            plot(30001:60:60000,t,'k')
        end
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(names{cnt},'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2savePics,names{cnt},'.png']);
    end
end



%% 20121023
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121023'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

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
% for i=1:39
%     subplot(5,8,i)
%     plot(mean(LinearFilter(:,1,49:2:60,i),3))
% end
onOff=[1 -1 1 1 -1 1 1 -1 1 1 -1 1 1 1 1 1 -1 1 1 1 1 1 1 1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 1 1];

path2savePics=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/FiringRateAndLinearFilter/'];
if ~exist(path2savePics,'dir')
    mkdir(path2savePics);
end
figure
set(gcf,'position',[560          17        1677         931])
nds='765432176543217654321'
for cnt=1:size(FiringRate,3)
    if onOff(cnt)~=0
        cc=1;
        for pt=0:2:5
             for i=(25+pt):24:(24*8+pt)
                subplot(3,7,cc)
                hold off
                if onOff(cnt)>0
                    col='r';
                else
                    col='b';
                end
                a=zeros(120000,4);
                p=1;
                for m=i:6:i+23
                    tmp=FiringRate(correctedProtocols(1,1,m):correctedProtocols(3600,1,m),m,cnt);
                    a(1:60000,p)=tmp(1:60000);
                    tmp=FiringRate(correctedProtocols(1,1,m+1):correctedProtocols(3600,1,m+1),m+1,cnt);
                    a(60001:120000,p)=tmp(1:60000);
                    p=p+1;
                end
                a=mean(a,2);
                plot(a,col)
                hold on
                
                t=mean(LinearFilter(:,1,i:6:(i+23),cnt),3);
                t=t-mean(t);
                t=t./sum(abs(t));
                ko=std(a(1:60000))/std(t);
                t=t*ko;
                tm=max(t);
                t=t-tm;
                plot(1:120:60000,t,col)
                
                t=mean(LinearFilter(:,1,(i+1):6:(i+23),cnt),3);
                t=t-mean(t);
                t=t./sum(abs(t));
                t=t*ko;
                t=t-tm;
                plot(60001:120:120000,t,'k')
                title(['ND',nds(cc)])
                cc=cc+1;
                axis tight
                axislimits(cc,1:2)=get(gca,'YLim');
            end
        end
        for cc=1:21
            subplot(3,7,cc)
            axis([0 120000 min(axislimits(:,1)) max(axislimits(:,2))])
            line([0,120000],[0,0],'color','k')
        end
        drawnow
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(names{cnt},'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2savePics,names{cnt},'.png']);
    end
end


%% 20121002
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121002'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

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
for i=1:30
    subplot(5,6,i)
    plot(mean(LinearFilter(:,1,49:2:60,i),3))
end
onOff=[-1 1 -1 -1 1 -1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 1 1 1 -1 -1 -1 1 1 1 -1 1];

path2savePics=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/FiringRateAndLinearFilter/'];
if ~exist(path2savePics,'dir')
    mkdir(path2savePics);
end
figure
set(gcf,'position',[560          17        1677         931])
nds='765432176543217654321'
for cnt=1:size(FiringRate,3)
    if onOff(cnt)~=0
        cc=1;
        for pt=0:2:5
             for i=(25+pt):24:(24*8+pt)
                subplot(3,7,cc)
                hold off
                if onOff(cnt)>0
                    col='r';
                else
                    col='b';
                end
                a=zeros(120000,4);
                p=1;
                for m=i:6:i+23
                    tmp=FiringRate(correctedProtocols(1,1,m):correctedProtocols(3600,1,m),m,cnt);
                    a(1:60000,p)=tmp(1:60000);
                    tmp=FiringRate(correctedProtocols(1,1,m+1):correctedProtocols(3600,1,m+1),m+1,cnt);
                    a(60001:120000,p)=tmp(1:60000);
                    p=p+1;
                end
                a=nanmean(a,2);
                plot(a,col)
                hold on
                
                t=nanmean(LinearFilter(:,1,i:6:(i+23),cnt),3);
                t=t-nanmean(t);
                t=t./sum(abs(t));
                if length(unique(t))>2
                    ko=std(a(1:60000))/std(t);
                    t=t*ko;
                    tm=max(t);
                    t=t-tm;
                    plot(1:120:60000,t,col)                    
                else
                    fl=1;
                end
                
                t=nanmean(LinearFilter(:,1,(i+1):6:(i+23),cnt),3);
                t=t-nanmean(t);
                t=t./sum(abs(t));
                if length(unique(t))>2
                    if fl
                        ko=std(a(60001:120000))/std(t);  
                        tm=max(t);
                    end
                    t=t*ko;
                    t=t-tm;
                    plot(60001:120:120000,t,'k')
                end
                fl=0;
                title(['ND',nds(cc)])
                cc=cc+1;
                axis tight
                axislimits(cc,1:2)=get(gca,'YLim');
            end
        end
        for cc=1:21
            subplot(3,7,cc)
            axis([0 120000 min(axislimits(:,1)) max(axislimits(:,2))])
            line([0,120000],[0,0],'color','k')
        end
        drawnow
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(names{cnt},'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2savePics,names{cnt},'.png']);
    end
end


%% 20120928
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120928'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

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
for i=1:24
    subplot(4,6,i)
    plot(mean(LinearFilter(:,1,25:2:48,i),3))
end
onOff=[-1 -1 1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 -1 1 -1 -1 -1 1 -1 -1 1];

path2savePics=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/FiringRateAndLinearFilter/'];
if ~exist(path2savePics,'dir')
    mkdir(path2savePics);
end
figure
set(gcf,'position',[560          17        1677         931])
nds='765432176543217654321'
fl=0;
for cnt=1:size(FiringRate,3)
    if onOff(cnt)~=0
        cc=1;
        for pt=0:2:5
             for i=(1+pt):24:(24*7+pt)
                subplot(3,7,cc)
                hold off
                if onOff(cnt)>0
                    col='r';
                else
                    col='b';
                end
                a=zeros(120000,4);
                p=1;
                for m=i:6:i+23
                    tmp=FiringRate(correctedProtocols(1,1,m):correctedProtocols(3600,1,m),m,cnt);
                    a(1:60000,p)=tmp(1:60000);
                    tmp=FiringRate(correctedProtocols(1,1,m+1):correctedProtocols(3600,1,m+1),m+1,cnt);
                    a(60001:120000,p)=tmp(1:60000);
                    p=p+1;
                end
                a=nanmean(a,2);
                plot(a,col)
                hold on
                
                t=nanmean(LinearFilter(:,1,i:6:(i+23),cnt),3);
                t=t-nanmean(t);
                t=t./sum(abs(t));
                if length(unique(t))>2
                    ko=std(a(1:60000))/std(t);
                    t=t*ko;
                    tm=max(t);
                    t=t-tm;
                    plot(1:120:60000,t,col)                    
                else
                    fl=1;
                end
                
                t=nanmean(LinearFilter(:,1,(i+1):6:(i+23),cnt),3);
                t=t-nanmean(t);
                t=t./sum(abs(t));
                if length(unique(t))>2
                    if fl
                        ko=std(a(60001:120000))/std(t);  
                        tm=max(t);
                    end
                    t=t*ko;
                    t=t-tm;
                    plot(60001:120:120000,t,'k')
                end
                fl=0;
                title(['ND',nds(cc)])
                cc=cc+1;
                axis tight
                axislimits(cc,1:2)=get(gca,'YLim');
            end
        end
        for cc=1:21
            subplot(3,7,cc)
            axis([0 120000 min(axislimits(:,1)) max(axislimits(:,2))])
            line([0,120000],[0,0],'color','k')
        end
        drawnow
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(names{cnt},'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2savePics,names{cnt},'.png']);
    end
end


%% 20130301
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130301'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

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
for i=1:33
    subplot(6,6,i)
    plot(mean(LinearFilter(:,1:2:end,5,i),2))
end
onOff=[-1 1 -1 1 1 -1 -1 -1 1 1 1 1 1 0 0 1 1 -1 1 -1 1 1 -1 -1 -1 1 -1 -1 -1 -1 1 1 1];

path2savePics=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/FiringRateAndLinearFilter/'];
if ~exist(path2savePics,'dir')
    mkdir(path2savePics);
end
figure
set(gcf,'position',[560          17        1677         931])
nds='87654321'
fl=0;
for cnt=1:size(FiringRate,3)
    if onOff(cnt)~=0
        cc=1;
        for i=1:2:16
            subplot(2,4,cc)
            if onOff(cnt)>0
                col='r';
            else
                col='b';
            end
            tmp=nanmean(FiringRate(correctedProtocols(1,1,it):correctedProtocols(3600,1,i),i:i+1,cnt),2);            
            tmp=tmp(1:60000);
            plot(tmp,col)
            hold on
            
            t=nanmean(LinearFilter(:,1:2:end,i:i+1,cnt),3);
            t=nanmean(t,2);
            t=t-nanmean(t);
            t=t./sum(abs(t));
            if length(unique(t))>2
                ko=std(tmp(1:10000))/std(t);
                t=t*ko;
                tm=max(t);
                t=t-tm;
                plot(1:60:30000,t,col)
            else
                fl=1;
            end
            
            t=nanmean(LinearFilter(:,2:2:end,i:i+1,cnt),3);
            t=nanmean(t,2);
            t=t-nanmean(t);
            t=t./sum(abs(t));
            if length(unique(t))>2
                if fl
                    ko=std(tmp(10001:20000))/std(t);
                    tm=max(t);
                end
                t=t*ko;
                t=t-tm;
                plot(30001:60:60000,t,'k')
            end
            fl=0;
            title(['ND',nds(cc)])
            axis tight
            axislimits(cc,1:2)=get(gca,'YLim');
            for mc=1:5
                lp=correctedProtocols(600*mc,1,i)-correctedProtocols(1,1,i);
                line([lp,lp],[0 300],'color',[0 0.7 0])
            end
            cc=cc+1;
        end
        for cc=1:8
            subplot(2,4,cc)
            axis([0 60000 min(axislimits(:,1)) max(axislimits(:,2))])
            line([0,60000],[0,0],'color','k')
            hold off
        end
        drawnow
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(names{cnt},'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2savePics,names{cnt},'.png']);
    end
end


%% 20130301_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130301_1'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

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
for i=1:34
    subplot(6,6,i)
    plot(mean(LinearFilter(:,1:2:end,5,i),2))
end
onOff=[1 -1 0 1 -1 -1 0 0 -1 -1 0 -1 0 -1 1 1 0 -1 -1 -1 1 -1 1 1 -1 1 0 1 -1 -1 1 -1 -1 0];

path2savePics=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/FiringRateAndLinearFilter/'];
if ~exist(path2savePics,'dir')
    mkdir(path2savePics);
end
figure
set(gcf,'position',[560          17        1677         931])
nds='87654321'
fl=0;
for cnt=1:size(FiringRate,3)
    if onOff(cnt)~=0
        cc=1;
        for i=1:2:16
            subplot(2,4,cc)
            if onOff(cnt)>0
                col='r';
            else
                col='b';
            end
            tmp=nanmean(FiringRate(correctedProtocols(1,1,i):correctedProtocols(3600,1,i),i:i+1,cnt),2);            
            tmp=tmp(1:60000);
            plot(tmp,col)
            hold on
            
            t=nanmean(LinearFilter(:,1:2:end,i:i+1,cnt),3);
            t=nanmean(t,2);
            t=t-nanmean(t);
            t=t./sum(abs(t));
            if length(unique(t))>2
                ko=std(tmp(1:10000))/std(t);
                t=t*ko;
                tm=max(t);
                t=t-tm;
                plot(1:60:30000,t,col)
            else
                fl=1;
            end
            
            t=nanmean(LinearFilter(:,2:2:end,i:i+1,cnt),3);
            t=nanmean(t,2);
            t=t-nanmean(t);
            t=t./sum(abs(t));
            if length(unique(t))>2
                if fl
                    ko=std(tmp(10001:20000))/std(t);
                    tm=max(t);
                end
                t=t*ko;
                t=t-tm;
                plot(30001:60:60000,t,'k')
            end
            fl=0;
            title(['ND',nds(cc)])
            axis tight
            axislimits(cc,1:2)=get(gca,'YLim');
            for mc=1:5
                lp=correctedProtocols(600*mc,1,i)-correctedProtocols(1,1,i);
                line([lp,lp],[0 300],'color',[0 0.7 0])
            end
            cc=cc+1;
        end
        for cc=1:8
            subplot(2,4,cc)
            axis([0 60000 min(axislimits(:,1)) max(axislimits(:,2))])
            line([0,60000],[0,0],'color','k')
            hold off
        end
        drawnow
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(names{cnt},'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2savePics,names{cnt},'.png']);
    end
end

%% 20130301_2
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130301_2'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

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
for i=1:36
    subplot(6,6,i)
    plot(mean(LinearFilter(:,1:2:end,5,i),2))
end
onOff=[-1 1 -1 1 1 -1 -1 -1 1 1 1 1 1 -1 -1 -1 -1 1 1 -1 1 -1 1 1 1 0 0 1 -1 1 0 1 -1 -1 0 -1];

path2savePics=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/FiringRateAndLinearFilter/'];
if ~exist(path2savePics,'dir')
    mkdir(path2savePics);
end
figure
set(gcf,'position',[560          17        1677         931])
nds='87654321'
fl=0;
for cnt=1:size(FiringRate,3)
    if onOff(cnt)~=0
        cc=1;
        for i=1:2:16
            subplot(2,4,cc)
            if onOff(cnt)>0
                col='r';
            else
                col='b';
            end
            tmp=nanmean(FiringRate(correctedProtocols(1,1,i):correctedProtocols(3600,1,i),i:i+1,cnt),2);            
            tmp=tmp(1:60000);
            plot(tmp,col)
            hold on
            
            t=nanmean(LinearFilter(:,1:2:end,i:i+1,cnt),3);
            t=nanmean(t,2);
            t=t-nanmean(t);
            t=t./sum(abs(t));
            if length(unique(t))>2
                ko=std(tmp(1:10000))/std(t);
                t=t*ko;
                tm=max(t);
                t=t-tm;
                plot(1:60:30000,t,col)
            else
                fl=1;
            end
            
            t=nanmean(LinearFilter(:,2:2:end,i:i+1,cnt),3);
            t=nanmean(t,2);
            t=t-nanmean(t);
            t=t./sum(abs(t));
            if length(unique(t))>2
                if fl
                    ko=std(tmp(10001:20000))/std(t);
                    tm=max(t);
                end
                t=t*ko;
                t=t-tm;
                plot(30001:60:60000,t,'k')
            end
            fl=0;
            title(['ND',nds(cc)])
            axis tight
            axislimits(cc,1:2)=get(gca,'YLim');
            for mc=1:5
                lp=correctedProtocols(600*mc,1,i)-correctedProtocols(1,1,i);
                line([lp,lp],[0 300],'color',[0 0.7 0])
            end
            cc=cc+1;
        end
        for cc=1:8
            subplot(2,4,cc)
            axis([0 60000 min(axislimits(:,1)) max(axislimits(:,2))])
            line([0,60000],[0,0],'color','k')
            hold off
        end
        drawnow
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(names{cnt},'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2savePics,names{cnt},'.png']);
    end
end


%% 20130302
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130302'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

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
for i=1:39
    subplot(5,8,i)
    plot(mean(LinearFilter(:,1:2:end,7,i),2))
end
onOff=[-1 1 0 0 0 0 0 -1 0 0 1 1 -1 1 1 1 -1 1 0 1 -1 1 1 -1 1 1 -1 1 -1 1 -1 1 1 -1 0 -1 1 -1 1];

path2savePics=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/FiringRateAndLinearFilter/'];
if ~exist(path2savePics,'dir')
    mkdir(path2savePics);
end
figure
set(gcf,'position',[560          17        1677         931])
nds='87654321'
fl=0;
for cnt=1:size(FiringRate,3)
    if onOff(cnt)~=0
        cc=1;
        for i=1:3:24
            subplot(2,4,cc)
            if onOff(cnt)>0
                col='r';
            else
                col='b';
            end
            tmp=nanmean(FiringRate(correctedProtocols(1,1,i):correctedProtocols(3600,1,i),i:i+2,cnt),2);            
            tmp=tmp(1:60000);
            plot(tmp,col)
            hold on
            
            t=nanmean(LinearFilter(:,1:2:end,i:i+2,cnt),3);
            t=nanmean(t,2);
            t=t-nanmean(t);
            t=t./sum(abs(t));
            if length(unique(t))>2
                ko=std(tmp(1:10000))/std(t);
                t=t*ko;
                tm=max(t);
                t=t-tm;
                plot(1:60:30000,t,col)
            else
                fl=1;
            end
            
            t=nanmean(LinearFilter(:,2:2:end,i:i+2,cnt),3);
            t=nanmean(t,2);
            t=t-nanmean(t);
            t=t./sum(abs(t));
            if length(unique(t))>2
                if fl
                    ko=std(tmp(10001:20000))/std(t);
                    tm=max(t);
                end
                t=t*ko;
                t=t-tm;
                plot(30001:60:60000,t,'k')
            end
            fl=0;
            title(['ND',nds(cc)])
            axis tight
            axislimits(cc,1:2)=get(gca,'YLim');
            for mc=1:5
                lp=correctedProtocols(600*mc,1,i)-correctedProtocols(1,1,i);
                line([lp,lp],[0 300],'color',[0 0.7 0])
            end
            cc=cc+1;
        end
        for cc=1:8
            subplot(2,4,cc)
            axis([0 60000 min(axislimits(:,1)) max(axislimits(:,2))])
            line([0,60000],[0,0],'color','k')
            hold off
        end
        drawnow
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(names{cnt},'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2savePics,names{cnt},'.png']);
    end
end

%% 20130302_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130302_1'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

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
for i=1:35
    subplot(6,6,i)
    plot(mean(LinearFilter(:,1:2:end,7,i),2))
end
onOff=[-1 0 0 -1 0 -1 1 -1 1 1 -1 -1 1 1 -1 1 1 1 -1 0 1 1 -1 0 1 -1 -1 1 -1 0 -1 -1 -1 0 -1];

path2savePics=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/FiringRateAndLinearFilter/'];
if ~exist(path2savePics,'dir')
    mkdir(path2savePics);
end
figure
set(gcf,'position',[560          17        1677         931])
nds='87654321'
fl=0;
for cnt=1:size(FiringRate,3)
    if onOff(cnt)~=0
        cc=1;
        for i=1:3:24
            subplot(2,4,cc)
            if onOff(cnt)>0
                col='r';
            else
                col='b';
            end
            tmp=nanmean(FiringRate(correctedProtocols(1,1,i):correctedProtocols(3600,1,i),i:i+2,cnt),2);            
            tmp=tmp(1:60000);
            plot(tmp,col)
            hold on
            
            t=nanmean(LinearFilter(:,1:2:end,i:i+2,cnt),3);
            t=nanmean(t,2);
            t=t-nanmean(t);
            t=t./sum(abs(t));
            if length(unique(t))>2
                ko=std(tmp(1:10000))/std(t);
                t=t*ko;
                tm=max(t);
                t=t-tm;
                plot(1:60:30000,t,col)
            else
                fl=1;
            end
            
            t=nanmean(LinearFilter(:,2:2:end,i:i+2,cnt),3);
            t=nanmean(t,2);
            t=t-nanmean(t);
            t=t./sum(abs(t));
            if length(unique(t))>2
                if fl
                    ko=std(tmp(10001:20000))/std(t);
                    tm=max(t);
                end
                t=t*ko;
                t=t-tm;
                plot(30001:60:60000,t,'k')
            end
            fl=0;
            title(['ND',nds(cc)])
            axis tight
            axislimits(cc,1:2)=get(gca,'YLim');
            for mc=1:5
                lp=correctedProtocols(600*mc,1,i)-correctedProtocols(1,1,i);
                line([lp,lp],[0 300],'color',[0 0.7 0])
            end
            cc=cc+1;
        end
        for cc=1:8
            subplot(2,4,cc)
            axis([0 60000 min(axislimits(:,1)) max(axislimits(:,2))])
            line([0,60000],[0,0],'color','k')
            hold off
        end
        drawnow
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title(names{cnt},'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2savePics,names{cnt},'.png']);
    end
end


%% 20130220
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130220'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

[k,m]=opt_subplots(length(units));
for i=1:length(units)    
    subplot(k,m,i)
    plot(mean(LinearFilter(:,1,29:37,i),3))
end
onOff=[-1 -1 -1 -1 -1 -1 1 1 0 1 -1 -1 1 1 1 1 0 -1 1 1 1 1 1 1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 1 -1 -1 -1 1 0 1 -1];

for i=1:length(units)    
    subplot(k,m,i)
    if onOff(i)>0
        col='r';
    else
        col='b';
    end
    plot(mean(LinearFilter(:,1,29:37,i),3),col)
end

%% 20130220_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130220_1'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

[k,m]=opt_subplots(length(units));
for i=1:length(units)    
    subplot(k,m,i)
    plot(mean(LinearFilter(:,1,28:36,i),3))
end
onOff=[1 -1 -1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 1 1 1 1];

for i=1:length(units)    
    subplot(k,m,i)
    if onOff(i)>0
        col='r';
    else
        col='b';
    end
    plot(mean(LinearFilter(:,1,29:37,i),3),col)
end


%% 20130224
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130224'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

[k,m]=opt_subplots(length(units));
for i=1:length(units)    
    subplot(k,m,i)
    plot(mean(LinearFilter(:,1,28:36,i),3))
end
onOff=[-1 1 1 1 -1 -1 1 1 1 -1 0 1 -1 -1 1 1 1 0 -1 -1 1 1 -1 1 -1 0 1 1 -1 1 1 -1 -1 -1 -1 -1 1 1 1 -1 -1 1 1];

for i=1:length(units)    
    subplot(k,m,i)
    if onOff(i)>0
        col='r';
    else
        col='b';
    end
    plot(mean(LinearFilter(:,1,29:37,i),3),col)
end


%% 20130225
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130225'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

[k,m]=opt_subplots(length(units));
for i=1:length(units)    
    subplot(k,m,i)
    plot(mean(LinearFilter(:,1,28:36,i),3))
end
onOff=[-1 1 1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 -1 -1 1 1 1 1 -1 1 -1 1 1 1 -1 -1 1 -1 1 1 -1 1 -1 1 0 1 1 -1 -1 1 1 1 1 1 1];

for i=1:length(units)
    subplot(k,m,i)
    if onOff(i)>0
        col='r';
    else
        col='b';
    end
    plot(mean(LinearFilter(:,1,29:37,i),3),col)
end


%% 20130226
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130226'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

[k,m]=opt_subplots(length(units));
for i=1:length(units)    
    subplot(k,m,i)
    plot(mean(LinearFilter(:,1,28:36,i),3))
end
onOff=[-1 -1 1 -1 -1 1 -1 -1 -1 0 -1 1 1 -1 1 -1 1 -1 1 -1 -1 -1 1 1 -1];

for i=1:length(units)    
    subplot(k,m,i)
    if onOff(i)>0
        col='r';
    else
        col='b';
    end
    plot(mean(LinearFilter(:,1,28:36,i),3),col)
end

%% 20130227
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130227'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

[k,m]=opt_subplots(length(units));
for i=1:length(units)    
    subplot(k,m,i)
    plot(mean(LinearFilter(:,1,13:16,i),3))
end
onOff=[1 1 0 -1 1 1 -1 0 0 0 1 1 -1 -1 0 1 -1 1 0 -1 0 0 0 0 1 1 0 0 -1 -1];

for i=1:length(units)    
    subplot(k,m,i)
    if onOff(i)>0
        col='r';
    elseif onOff(i)<0
        col='b';
    else
        col='g';
    end
    plot(mean(LinearFilter(:,1,21:24,i),3),col)
end

for i=1:36
    subplot(6,6,i)
    plot(LinearFilter(:,1,i,27))
end


%% 20130227
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130227'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

[k,m]=opt_subplots(length(units));
for i=1:length(units)    
    subplot(k,m,i)
    plot(mean(LinearFilter(:,1,13:16,i),3))
end
onOff=[1 1 0 -1 1 1 -1 0 0 0 1 1 -1 -1 0 1 -1 1 0 -1 0 0 0 0 1 1 0 0 -1 -1];

for i=1:length(units)    
    subplot(k,m,i)
    if onOff(i)>0
        col='r';
    elseif onOff(i)<0
        col='b';
    else
        col='g';
    end
    plot(mean(LinearFilter(:,1,21:24,i),3),col)
end

for i=1:36
    subplot(6,6,i)
    plot(LinearFilter(:,1,i,27))
end

%% 20120329
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120329'
codeWord='ffflicker'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

[k,m]=opt_subplots(length(units));
for i=1:length(units)    
    subplot(k,m,i)
    plot(mean(LinearFilter(:,1,9:13,i),3))
end
onOff=[1 1 -1 1 -1 1 -1 -1 1 -1 0 1 1 0 -1 1 0 0 0 1 -1 1 1 1 -1 0 0 1 0 0 1 1 0 1 -1 -1 1 -1 -1 -1 0 1 1 1 1 1 1 -1 1 1 -1 -1];

for i=1:length(units)    
    subplot(k,m,i)
    if onOff(i)>0
        col='r';
    elseif onOff(i)<0
        col='b';
    else
        col='g';
    end
    plot(mean(LinearFilter(:,1,9:13,i),3),col)
end

%% 20120627
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120627'
codeWord='ffflicker'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

[k,m]=opt_subplots(length(units));
for i=1:length(units)    
    subplot(k,m,i)
    plot(mean(LinearFilter(:,1,13:18,i),3))
end
onOff=[0 -1 1 -1 1 0 1 -1 -1 1 1 1 1 -1 1 1 -1 0];

for i=1:length(units)    
    subplot(k,m,i)
    if onOff(i)>0
        col='r';
    elseif onOff(i)<0
        col='b';
    else
        col='g';
    end
    plot(mean(LinearFilter(:,1,13:18,i),3),col)
end

%% 20120714
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120714'
codeWord='ffflicker'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')

[k,m]=opt_subplots(length(units));
for i=1:length(units)    
    subplot(k,m,i)
    plot(mean(LinearFilter(:,1,7:12,i),3))
end
onOff=[0 1 1 1 0 0 1 1 0 0 0 1 0 1 -1];

for i=1:length(units)    
    subplot(k,m,i)
    if onOff(i)>0
        col='r';
    elseif onOff(i)<0
        col='b';
    else
        col='g';
    end
    plot(mean(LinearFilter(:,1,7:12,i),3),col)
end
