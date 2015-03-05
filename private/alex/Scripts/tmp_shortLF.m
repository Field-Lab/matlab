%-- 2/21/2013 1:55 PM --%
clear
date='20130226';
unitspath=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
units=dir([unitspath,'*FFFlicker.mat']);
path2save=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
if ~exist(path2save,'file')
    mkdir(path2save);
end
load([unitspath,units(1).name])
filter_length=500; %filter length in ms
lseq=[];hseq=[];%m=[];
for i=1:length(spike_info.name_info)
    if ~isempty(cell2mat(regexp(spike_info.name_info(i),'LCseq')))
        lseq=[lseq i];
    else
        %         a=cell2mat(regexp(spike_info.name_info(i),'ND'))+2;
        %         if length(a)>1
        %             k=str2num(spike_info.name_info{i}(a(1)))+str2num(spike_info.name_info{i}(a(2)))
        %         else
        %             k=str2num(spike_info.name_info{i}(a));
        %         end
        %         m=[m k];
        hseq=[hseq i];
    end
end
for units_cnt=1:length(units)
    units(units_cnt).name
    if ~exist([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter.mat'],'file')
        load([unitspath,units(units_cnt).name])
        trials=size(spike_info.flip_times,1);
        HighSpikeCount=zeros(1,length(hseq));
        LowSpikeCount=zeros(1,length(lseq));
        HighFilters=zeros(filter_length,length(hseq));
        LowFilters=zeros(filter_length,length(lseq));
        cntL=1;
        cntH=1;
        for trial_counter=1:trials
            flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
            pxls=(((spike_info.flip_times{trial_counter}(:,2))-30)/30); % pxs in ms, normalized -1:1
            spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
            spikes(spikes<filter_length|spikes>flips(end))=[];
            if isempty(spikes)
                if ~isempty(find(hseq==trial_counter, 1))
                    cntH=cntH+1;
                else
                    cntL=cntL+1;
                end
            else
                %         if trial_counter==146 % for 20120920
                %             flips=spike_info.flip_times{trial_counter-1}(:,1); % flips in ms, rounded
                %         end
                tmp=zeros(flips(end)+1,2); % 1 column - color, 2 column - spikes
                tmp([1; flips(1:end-1)],1)=pxls;
                % fill stimulus down to ms with brightness values where spikes
                % happen
                for i=1:50
                    subscr=(flips(1:end-1)+i);
                    subscr(tmp(flips(1:end-1)+i,1)~=0)=[];
                    tmp(subscr,1)=tmp(subscr-1,1);
                end
                tmp(spikes,2)=1;
                n=zeros(length(spikes),filter_length); % filter
                for i=1:filter_length
                    n(:,i)=tmp(spikes-i+1,1);
                    %     n(i)=sum(tmp(c-i+1,1))/all_spikes;
                end
                n=[spikes, n];
                % plot(b)
                tmp=n(n(:,1)<flips(3600)&n(:,1)>filter_length,2:end);
                if ~isempty(find(hseq==trial_counter, 1))
                    HighFilters(1:filter_length,cntH)=sum(tmp);
                    HighSpikeCount(cntH)=size(tmp,1);
                    cntH=cntH+1;
                else
                    LowFilters(1:filter_length,cntL)=sum(tmp);
                    LowSpikeCount(cntL)=size(tmp,1);
                    cntL=cntL+1;
                end
            end
        end
        save([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter.mat'],'HighFilters','HighSpikeCount','LowFilters','LowSpikeCount')
    end
end

% plotting special for 20130224
cd('S:\user\alexandra\scripts')
clear
date='20130224';
filter_length=500;
nd='987654321';
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
path2save=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF_plots\'];
if ~exist(path2save,'file')
    mkdir(path2save);
end
figure(1)
set(gcf,'Position', [1 31 1680 946])
[rows,cols]=opt_subplots(length(nd));
filters_cnt=1
load([filterspath,filters(filters_cnt).name])
minn=min(min(HighFilters(:,1:2:end)));
maxx=max(max(HighFilters(:,1:2:end)));
cnt=1;
for i=[1,10:9:82]
    k=subplot(rows,cols,cnt);
    plot(HighFilters(:,i:i+8))
    line([0 500],[0 0],'color',[1 1 1]*0.5)
    hleg=legend(int2str(HighSpikeCount(i:i+8)'));
    legend('boxoff')
    set(hleg,'Fontsize',8)
    title(['ND',nd(cnt)])
    axis tight
    %         axis([0 500 -100 150])
    %axis([0 500 minn maxx])
    cnt=cnt+1;
end
k
i
[1,10:9:81]
256/4
minn=min(min(HighFilters(:,1:2:end)));
maxx=max(max(HighFilters(:,1:2:end)));
cnt=1;
for i=[1,10:9:81]
k=subplot(rows,cols,cnt);
plot(HighFilters(:,i:i+8))
line([0 500],[0 0],'color',[1 1 1]*0.5)
hleg=legend(int2str(HighSpikeCount(i:i+8)'));
legend('boxoff')
set(hleg,'Fontsize',8)
title(['ND',nd(cnt)])
axis tight
%         axis([0 500 -100 150])
%axis([0 500 minn maxx])
cnt=cnt+1;
end
[1,11:9:81]
minn=min(min(HighFilters(:,1:2:end)));
maxx=max(max(HighFilters(:,1:2:end)));
cnt=1;
for i=[1,11:9:81]
k=subplot(rows,cols,cnt);
plot(HighFilters(:,i:i+8))
line([0 500],[0 0],'color',[1 1 1]*0.5)
hleg=legend(int2str(HighSpikeCount(i:i+8)'));
legend('boxoff')
set(hleg,'Fontsize',8)
title(['ND',nd(cnt)])
axis tight
%         axis([0 500 -100 150])
%axis([0 500 minn maxx])
cnt=cnt+1;
end
title([filters(filters_cnt).name(1:end-4),' ',int2str(filter_length),'ms   High CONTRAST'],'FontSize',15,'FontWeight','bold','Interpreter','none')
for filters_cnt=1:length(filters)
if ~exist([path2save,filters(filters_cnt).name(1:end-3),'emf'],'file')
load([filterspath,filters(filters_cnt).name])
minn=min(min(HighFilters(:,1:2:end)));
maxx=max(max(HighFilters(:,1:2:end)));
cnt=1;
for i=[1,11:9:81]
k=subplot(rows,cols,cnt);
plot(HighFilters(:,i:i+8))
line([0 500],[0 0],'color',[1 1 1]*0.5)
hleg=legend(int2str(HighSpikeCount(i:i+8)'));
legend('boxoff')
set(hleg,'Fontsize',8)
title(['ND',nd(cnt)])
axis tight
%         axis([0 500 -100 150])
%axis([0 500 minn maxx])
cnt=cnt+1;
end
subplot('Position',[0.5 0.95 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title([filters(filters_cnt).name(1:end-4),' ',int2str(filter_length),'ms   High CONTRAST'],'FontSize',15,'FontWeight','bold','Interpreter','none')
saveas(gcf,[path2save,filters(filters_cnt).name(1:end-3),'emf'])
end
end
figure(1)
set(gcf,'Position', [1 31 1680 946])
[rows,cols]=opt_subplots(length(nd));
% simple filter plot
for filters_cnt=1:length(filters)
if ~exist([path2save,filters(filters_cnt).name(1:end-3),'emf'],'file')
load([filterspath,filters(filters_cnt).name])
minn=min(min(HighFilters(:,1:2:end)));
maxx=max(max(HighFilters(:,1:2:end)));
cnt=1;
for i=[1,11:9:81]
    k=subplot(rows,cols,cnt);
    plot(HighFilters(:,i:i+8))
    line([0 500],[0 0],'color',[1 1 1]*0.5)
    hleg=legend(int2str(HighSpikeCount(i:i+8)'));
    legend('boxoff')
    set(hleg,'Fontsize',8)
    title(['ND',nd(cnt)])
    axis tight
    %         axis([0 500 -100 150])
    %axis([0 500 minn maxx])
    cnt=cnt+1;
end
subplot('Position',[0.5 0.95 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title([filters(filters_cnt).name(1:end-4),' ',int2str(filter_length),'ms   High CONTRAST'],'FontSize',15,'FontWeight','bold','Interpreter','none')
saveas(gcf,[path2save,filters(filters_cnt).name(1:end-3),'emf'])
end
end
clear
path2data=['S:\user\alexandra\MEA_data\',date,'\'];
unitsPath=dir([path2data,'easy_formatted_units\*FFFlash*']);
path2save=[path2data,'Full_info\'];
if ~exist(path2save,'file')
mkdir(path2save);
end
load([path2data,'easy_formatted_units\',unitsPath(1).name])
trials=size(spike_info.name_info,1);
unitsPath
clear
date='20130220'
path2data=['S:\user\alexandra\MEA_data\',date,'\'];
path2data
unitsPath=dir([path2data,'easy_formatted_units\*FFFlash*']);
[path2data,'easy_formatted_units\*FFFlash*']
path2data=['S:\data\alexandra\MEA_data\',date,'\'];
unitsPath=dir([path2data,'easy_formatted_units\*FFFlash*']);
path2save=[path2data,'FFFlash\'];
if ~exist(path2save,'file')
mkdir(path2save);
end
load([path2data,'easy_formatted_units\',unitsPath(1).name])
trials=size(spike_info.name_info,1);
trials
nds='987654321';
filterspath=['F:\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
unit=1
load([path2fit,fits(unit).name]);
load([path2data,'easy_formatted_units\',unitsPath(unit).name,'___spike_info_FFFlash'])
name
unitsPath(unit).name
load([path2data,'easy_formatted_units\',unitsPath(unit).name])
white_flash=zeros(trials,4501);
black_flash=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash(trial,1:4501)=white_flash(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash(trial,1:4501)=black_flash(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
i=1
subplot(3,3,i)
hold on
area(1:500,white_flash(i+1,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(i+1,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(i+1,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(i+1,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(i+1,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(i+1,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(i+1,1:4500) black_flash(i+1,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
for i=1:9
subplot(3,3,i)
hold on
area(1:500,white_flash(i+1,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(i+1,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(i+1,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(i+1,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(i+1,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(i+1,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(i+1,1:4500) black_flash(i+1,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
end
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
for i=1:9
subplot(3,3,i)
hold on
area(1:500,white_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(i,1:4500) black_flash(i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
end
trials
white_flash=zeros(trials-1,4501);
black_flash=zeros(trials-1,4501);
trial=1
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
white_flash=zeros(trials,4501);
black_flash=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash(trial,1:4501)=white_flash(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash(trial,1:4501)=black_flash(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
white_flash=white_flash/5;
black_flash=black_flash/5;
white_flash
plot(white_flash)
plot(white_flash#)
plot(white_flash')
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
for i=[1 11:9:81]
white_flash(i)=mean(white_flash_tmp(i:i+8,:));
black_flash(i)=mean(black_flash_tmp(i:i+8,:));
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
for i=1:9
subplot(3,3,i)
hold on
area(1:500,white_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(i,1:4500) black_flash(i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(unitsPath(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
[path2save,unitsPath(unit).name(1:end-20),'.emf']
[path2save,unitsPath(unit).name(1:end-25),'.emf']
load([path2data,'easy_formatted_units\',unitsPath(unit).name])
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
for i=1:9
subplot(3,3,i)
hold on
area(1:500,white_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(i,1:4500) black_flash(i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
title(['ND',int2str(10-i)])
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(unitsPath(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
saveas(gcf,[path2save,unitsPath(unit).name(1:end-25),'.emf']);
% FULL FIELD FLASH
load([path2data,'easy_formatted_units\',unitsPath(unit).name])
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
for i=1:9
subplot(9,3,i)
hold on
area(1:500,white_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(i,1:4500) black_flash(i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
title(['ND',int2str(10-i)])
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(unitsPath(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
load([path2data,'easy_formatted_units\',unitsPath(unit).name])
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=1;
for i=1:9
subplot(9,3,kk)
hold on
area(1:500,white_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(i,1:4500) black_flash(i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
title(['ND',int2str(10-i)])
kk=kk+3;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(unitsPath(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
kk=1;
for i=1:9
subplot(9,3,kk)
hold on
area(1:500,white_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(i,1:4500) black_flash(i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(10-i)])
kk=kk+3;
end
load([path2data,'easy_formatted_units\',unitsPath(unit).name])
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=1;
for i=1:9
subplot(9,3,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(10-i)])
kk=kk+3;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(unitsPath(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
kk=1;
for i=1:9
subplot(9,3,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+3;
end
load([filterspath,filters(unit).name])
filters
LF
LF=HighFilters([1:9 11:82],:);
i=1
9*(i-1)+1:9*i
i=2
9*(i-1)+1:9*i
i=1
9*(10-i-1)+1:9*(10-i)
i=9
9*(10-i-1)+1:9*(10-i)
kk=2;
for i=1:9
subplot(9,3,kk)
plot(LF(9*(10-i-1)+1:9*(10-i))')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+3;
end
LF=HighFilters(:,[1:9 11:82]);
LF=HighFilters(:,[1:9 11:82])';
kk=2;
for i=1:9
subplot(9,3,kk)
plot(LF(9*(10-i-1)+1:9*(10-i))')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+3;
end
kk=2;
for i=1:9
subplot(9,3,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+3;
end
path2moment=['S:\data\alexandra\MEA_data\',date,'\'];
momentPath=dir([path2moment,'easy_formatted_units\*Moment*']);
momentPath
load([path2moment,'easy_formatted_units\',momentPath(unit).name])
white_flash_tmp=zeros(trials,6000);
black_flash_tmp=zeros(trials,6000);
trial=1
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips
diff(flips(:,1))
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
st=1
flips(st)-
flips(st)
flips
flips(st):flips(st+4)
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
load([path2moment,'easy_formatted_units\',momentPath(unit).name])
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
trial
st
load([path2moment,'easy_formatted_units\',momentPath(unit).name])
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
maxx
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=2;
for i=1:9
subplot(9,3,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+3;
end
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
maxx
load([path2moment,'easy_formatted_units\',momentPath(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=2;
for i=1:9
subplot(9,3,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+3;
end
% FULL FIELD FLASH
load([path2data,'easy_formatted_units\',unitsPath(unit).name])
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=1;
for i=1:9
subplot(9,3,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+3;
end
% Moment
load([path2moment,'easy_formatted_units\',momentPath(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=2;
for i=1:9
subplot(9,3,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+3;
end
kk=3;
for i=1:9
subplot(9,3,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+3;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(unitsPath(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
clc
clear
date='20130220'
path2data=['S:\data\alexandra\MEA_data\',date,'\'];
unitsPath=dir([path2data,'easy_formatted_units\*FFFlash*']);
path2moment=['S:\data\alexandra\MEA_data\',date,'\'];
momentPath=dir([path2moment,'easy_formatted_units\*Moment*']);
path2save=[path2data,'full_info\'];
if ~exist(path2save,'file')
mkdir(path2save);
end
load([path2data,'easy_formatted_units\',unitsPath(1).name])
trials=size(spike_info.name_info,1);
nds='987654321';
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(fits)
% FULL FIELD FLASH
load([path2data,'easy_formatted_units\',unitsPath(unit).name])
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=1;
for i=1:9
subplot(9,3,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+3;
end
% Moment
load([path2moment,'easy_formatted_units\',momentPath(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=2;
for i=1:9
subplot(9,3,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+3;
end
kk=3;
for i=1:9
subplot(9,3,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+3;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(unitsPath(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
saveas(gcf,[path2save,unitsPath(unit).name(1:end-25),'.emf']);
end
for unit=1:length(unitsPath)
% FULL FIELD FLASH
load([path2data,'easy_formatted_units\',unitsPath(unit).name])
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=1;
for i=1:9
subplot(9,3,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+3;
end
% Moment
load([path2moment,'easy_formatted_units\',momentPath(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=2;
for i=1:9
subplot(9,3,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+3;
end
kk=3;
for i=1:9
subplot(9,3,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+3;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(unitsPath(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
saveas(gcf,[path2save,unitsPath(unit).name(1:end-25),'.emf']);
end
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\];
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
clear
date='20130220'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);
moments=dir([path2moment,'*Moment*']);
moments=dir([path2data,'*Moment*']);
nds=dir([path2data,'*nd*']);
path2save=['S:\data\alexandra\MEA_data\',date,'\full_info\'];
path2save
if ~exist(path2save,'file')
mkdir(path2save);
end
load([path2data,flashes(1).name])
trials=size(spike_info.name_info,1);
nds='987654321';
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
unit=1
load([path2data,flashes(unit).name])
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=1;
kk=2;
for i=1:9
subplot(9,4,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=3;
for i=1:9
subplot(9,4,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
kk=4;
for i=1:9
subplot(9,4,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
load([path2data,nds(unit).name])
unit
nds
nds=dir([path2data,'*nd*']);
load([path2data,nds(unit).name])
spike_info.spike_times
trial=[1
]
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
flips
spike_info.flip_times(trial)
cell2mat(spike_info.flip_times(trial))
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
flips(end)
spike_info.flip_times(trial)
spike_info.flip_times(trial){1}
spike_info.flip_times(trial)(1)
cell2mat(spike_info.flip_times(trial))
flips=flips(1,end);
conv=convolved(spikes,40,flips);
conv=conv(121:end-120);
nd_trials=[1 3:10];
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
flips=flips(1,end);
conv=convolved(spikes,40,flips);
conv=conv(121:end-120);
plot(conv)
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
conv
nd_trials
nd_trials(i)
cell2mat(spike_info.spike_times(nd_trials(i)))
flips
flips=flips(1,end);
flips
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
flips
flips=flips(end,1);
flips
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
flips=flips(end,1);
conv=convolved(spikes,40,flips);
conv=conv(121:end-120);
plot(conv)
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
nd_trials=[10:-1:3 1];
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
flips=flips(end,1);
conv=convolved(spikes,40,flips);
conv=conv(121:end-120);
plot(conv)
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
flips
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
flips
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
flips=flips(end,1);
conv=convolved(spikes,40,flips);
conv=conv(121:end-120);
plot(conv)
line([flips(1,2),flips(1,2)],[0,max(conv)],'color','k')
line([flips(6,2),flips(6,2)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv)
line([flips(1,2),flips(1,2)],[0,max(conv)],'color','k')
line([flips(6,2),flips(6,2)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
max(conv)
kk=1;
i=1
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv)
line([flips(1,2),flips(1,2)],[0,max(conv)],'color','k')
[flips(1,2),flips(1,2)]
flips
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv)
line([flips(1,1),flips(1,1)],[0,max(conv)],'color','k')
line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
flips(1,1)
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv)
line([flips(1,1)-1000,flips(1,1)],[0,max(conv)],'color','k')
line([flips(6,1)-1000,flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv)
line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
line([flips(6,1)-1000,flips(6,1)-1000],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv)
line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv,'Linewidth',2)
line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
for unit=1:length(unitsPath)
%nd
load([path2data,nds(unit).name])
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv,'Linewidth',2)
line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
for trial=[1 3:10]
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(1,end);
conv=convolved(spikes,40,flips);
conv=conv(121:end-120);
end
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=2;
for i=1:9
subplot(9,4,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% Moment
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=3;
for i=1:9
subplot(9,4,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% LF
kk=4;
for i=1:9
subplot(9,4,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(unitsPath(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
saveas(gcf,[path2save,unitsPath(unit).name(1:end-25),'.emf']);
end
for unit=1:length(unitsPath)
%nd
load([path2data,nds(unit).name])
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv,'Linewidth',2)
line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
for trial=[1 3:10]
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(1,end);
conv=convolved(spikes,40,flips);
conv=conv(121:end-120);
end
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=2;
for i=1:9
subplot(9,4,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% Moment
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=3;
for i=1:9
subplot(9,4,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% LF
kk=4;
for i=1:9
subplot(9,4,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(path2data(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
saveas(gcf,[path2save,path2data(unit).name(1:end-25),'.emf']);
end
for unit=1:length(path2data)
%nd
load([path2data,nds(unit).name])
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv,'Linewidth',2)
line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
for trial=[1 3:10]
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(1,end);
conv=convolved(spikes,40,flips);
conv=conv(121:end-120);
end
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=2;
for i=1:9
subplot(9,4,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% Moment
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=3;
for i=1:9
subplot(9,4,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% LF
kk=4;
for i=1:9
subplot(9,4,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(path2data(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
saveas(gcf,[path2save,path2data(unit).name(1:end-25),'.emf']);
end
unit
load([path2data,nds(unit).name])
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv,'Linewidth',2)
line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
for unit=1:length(path2data)
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
%nd
load([path2data,nds(unit).name])
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv,'Linewidth',2)
line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=2;
for i=1:9
subplot(9,4,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% Moment
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=3;
for i=1:9
subplot(9,4,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% LF
kk=4;
for i=1:9
subplot(9,4,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(path2data(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
saveas(gcf,[path2save,path2data(unit).name(1:end-25),'.emf']);
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(path2data(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
path2data(unit)
path2data
for unit=1:length(path2data)
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
%nd
load([path2data,nds(unit).name])
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(nd_trials(i)));
flips=cell2mat(spike_info.flip_times(nd_trials(i)));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv,'Linewidth',2)
line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=2;
for i=1:9
subplot(9,4,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% Moment
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=3;
for i=1:9
subplot(9,4,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% LF
kk=4;
for i=1:9
subplot(9,4,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(flashes(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
saveas(gcf,[path2save,flashes(unit).name(1:end-25),'.emf']);
end
length(path2data)
unit
%get linear filter
clear
date='20130220_1';
unitspath=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
units=dir([unitspath,'*FFFlicker.mat']);
path2save=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
if ~exist(path2save,'file')
mkdir(path2save);
end
load([unitspath,units(1).name])
filter_length=500; %filter length in ms
lseq=[];hseq=[];%m=[];
for i=1:length(spike_info.name_info)
if ~isempty(cell2mat(regexp(spike_info.name_info(i),'LCseq')))
lseq=[lseq i];
else
%         a=cell2mat(regexp(spike_info.name_info(i),'ND'))+2;
%         if length(a)>1
%             k=str2num(spike_info.name_info{i}(a(1)))+str2num(spike_info.name_info{i}(a(2)))
%         else
%             k=str2num(spike_info.name_info{i}(a));
%         end
%         m=[m k];
hseq=[hseq i];
end
end
for units_cnt=1:length(units)
units(units_cnt).name
if ~exist([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter.mat'],'file')
load([unitspath,units(units_cnt).name])
trials=size(spike_info.flip_times,1);
HighSpikeCount=zeros(1,length(hseq));
LowSpikeCount=zeros(1,length(lseq));
HighFilters=zeros(filter_length,length(hseq));
LowFilters=zeros(filter_length,length(lseq));
cntL=1;
cntH=1;
for trial_counter=1:trials
flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
pxls=(((spike_info.flip_times{trial_counter}(:,2))-30)/30); % pxs in ms, normalized -1:1
spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
spikes(spikes<filter_length|spikes>flips(end))=[];
if isempty(spikes)
if ~isempty(find(hseq==trial_counter, 1))
cntH=cntH+1;
else
cntL=cntL+1;
end
else
%         if trial_counter==146 % for 20120920
%             flips=spike_info.flip_times{trial_counter-1}(:,1); % flips in ms, rounded
%         end
tmp=zeros(flips(end)+1,2); % 1 column - color, 2 column - spikes
tmp([1; flips(1:end-1)],1)=pxls;
% fill stimulus down to ms with brightness values where spikes
% happen
for i=1:50
subscr=(flips(1:end-1)+i);
subscr(tmp(flips(1:end-1)+i,1)~=0)=[];
tmp(subscr,1)=tmp(subscr-1,1);
end
tmp(spikes,2)=1;
n=zeros(length(spikes),filter_length); % filter
for i=1:filter_length
n(:,i)=tmp(spikes-i+1,1);
%     n(i)=sum(tmp(c-i+1,1))/all_spikes;
end
n=[spikes, n];
% plot(b)
tmp=n(n(:,1)<flips(3600)&n(:,1)>filter_length,2:end);
if ~isempty(find(hseq==trial_counter, 1))
HighFilters(1:filter_length,cntH)=sum(tmp);
HighSpikeCount(cntH)=size(tmp,1);
cntH=cntH+1;
else
LowFilters(1:filter_length,cntL)=sum(tmp);
LowSpikeCount(cntL)=size(tmp,1);
cntL=cntL+1;
end
end
end
save([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter.mat'],'HighFilters','HighSpikeCount','LowFilters','LowSpikeCount')
end
end
clear
date='20130220_1'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);
moments=dir([path2data,'*Moment*']);
nds=dir([path2data,'*nd*']);
path2save=['S:\data\alexandra\MEA_data\',date,'\full_info\'];
if ~exist(path2save,'file')
mkdir(path2save);
end
load([path2data,flashes(1).name])
trials=size(spike_info.name_info,1);
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(flashes)
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
%nd
load([path2data,nds(unit).name])
kk=1;
for i=1:9
subplot(9,4,kk)
spikes=cell2mat(spike_info.spike_times(10-i));
flips=cell2mat(spike_info.flip_times(10-i));
conv=convolved(spikes,40,flips(end,1));
conv=conv(121:end-120);
plot(conv,'Linewidth',2)
line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
title(['ND',int2str(i)])
axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
load([filterspath,filters(unit).name])
LF=HighFilters';
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=1:9:81
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=2;
for i=1:9
subplot(9,4,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% Moment
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=1:9:81
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
kk=3;
for i=1:9
subplot(9,4,kk)
hold on
plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
% LF
kk=4;
for i=1:9
subplot(9,4,kk)
plot(LF(9*(10-i-1)+1:9*(10-i),:)')
line([0 500],[0 0],'color',[1 1 1]*0.5)
title(['ND',int2str(i)])
axis tight
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
kk=kk+4;
end
subplot('Position',[0.5 0.97 0.00001 0.00001])
set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
title(flashes(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
saveas(gcf,[path2save,flashes(unit).name(1:end-25),'.emf']);
end
clear
date='20130220'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
unit=1
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
i=1
meanFilter=mean(LF(i:i+8,:));
filters(unit).name
meanFilter=zeros(9,500);
cnt=1;
for i=1:9:80
meanFilter(cnt,:)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
plot(meanFilter)
plot(meanFilter')
meanFilter=zeros(length(filters),500,9);
length(filters)
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
plot(meanFilter(:,:,1))
plot(meanFilter(:,:,1)')
for i=1:9
subplot(3,3,i)
plot(meanFilter(:,:,1)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(:,:,i)')
end
tmp=mean(LF(i:i+8,:));
tmp
max(abs(tmp))
meanFilter=zeros(length(filters),500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
tmp=mean(LF(i:i+8,:));
meanFilter(unit,:,cnt)=tmp/max(abs(tmp));
cnt=cnt+1;
end
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(:,:,i)')
end
cool
std(meanFilter(4,:,:),0,3)
meanFilter(4,:,:)
std(reshape(meanFilter(4,:,:),500,9),0,2)
std(reshape(meanFilter(4,:,:),500,9))
meanFilter(4,:,:)
std(meanFilter(:,:,4))
meanFilter(:,:,4)
std(meanFilter(:,:,4),0,2)
for i=1:9
subplot(3,3,i)
plot(meanFilter(:,:,i)')
end
snr=std(meanFilter(:,60:300,4),0,2)./std(meanFilter(:,1:20,4),0,2);
snr
plot(snr)
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>20,:,i)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>25,:,i)')
end
snr>25
sum(snr>25)
sum(snr>30)
sum(snr>20)
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>25,:,i)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>25,:,i)')
end
mean(meanFilter(:,130:140,4))
meanFilter=zeros(length(filters),500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
mean(meanFilter(:,130:140,4))
mean(meanFilter(:,130:140,4),2)
meanFilter(unit,:,:)
max(meanFilter(unit,:,:))
a=reshape(max(meanFilter(unit,:,:)),9,1)
meanFilter(unit,:,:)
reshape(max(meanFilter(unit,:,:)),9,1)
max(meanFilter(unit,:,:))
meanFilter(unit,:,:)
tmp=reshape(meanFilter(unit,:,:),500,9)';
plot(tmp)
plot(tmp')
tmp=reshape(meanFilter(unit,:,:),500,9)'./repmat(a,1,500);
plot(tmp')
mean(meanFilter(unit,130:140,4),2)>0
meanFilter(unit,:,:)
meanFilter(unit,:,:)=tmp';
meanFilter=zeros(length(filters),500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
for unit=1:length(filters)
if mean(meanFilter(unit,130:140,4),2)>0
a=reshape(max(meanFilter(unit,:,:)),9,1);
else
a=reshape(min(meanFilter(unit,:,:)),9,1);
end
tmp=reshape(meanFilter(unit,:,:),500,9)'./repmat(abs(a),1,500);
meanFilter(unit,:,:)=tmp';
end
snr=std(meanFilter(:,60:300,4),0,2)./std(meanFilter(:,1:20,4),0,2);
plot(snr)
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>25,:,i)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30,:,i)')
end
mean(meanFilter(unit,130:140,4),2)
meanFilter(:,130:140,4)
mean(meanFilter(:,130:140,4),2)
onOff=mean(meanFilter(:,130:140,4),2)
onOff=mean(meanFilter(:,130:140,4),2);
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30&onOff>0,:,i)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30&onOff<0,:,i)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30&onOff<0,1:300,i)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30&onOff>0,1:300,i)')
end
clear
date='20130220'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(length(filters),500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters';
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
clear
date='20130220'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(length(filters),500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
clear
date='20130220'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(length(filters),500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
unit=1
load([filterspath,filters(unit).name])
LF=HighFilters';
date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters';
for i=1:9:80
meanFilter(size(meanFilter,1)+1,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
clear
date='20130220'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(46+28,500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters';
for i=1:9:80
meanFilter(unit+46,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
for unit=1:size(meanFilter,1)
if mean(meanFilter(unit,130:140,4),2)>0
a=reshape(max(meanFilter(unit,:,:)),9,1);
else
a=reshape(min(meanFilter(unit,:,:)),9,1);
end
tmp=reshape(meanFilter(unit,:,:),500,9)'./repmat(abs(a),1,500);
meanFilter(unit,:,:)=tmp';
end
unit=1
mean(meanFilter(unit,130:140,4),2)>0
a=reshape(min(meanFilter(unit,:,:)),9,1);
meanFilter(unit,:,:)
unit
min(meanFilter(unit,:,:))
clear
date='20130220'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(46+28,500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters';
cnt=1;
for i=1:9:80
meanFilter(unit+46,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
for unit=1:size(meanFilter,1)
if mean(meanFilter(unit,130:140,4),2)>0
a=reshape(max(meanFilter(unit,:,:)),9,1);
else
a=reshape(min(meanFilter(unit,:,:)),9,1);
end
tmp=reshape(meanFilter(unit,:,:),500,9)'./repmat(abs(a),1,500);
meanFilter(unit,:,:)=tmp';
end
snr=std(meanFilter(:,60:300,4),0,2)./std(meanFilter(:,1:20,4),0,2);
onOff=mean(meanFilter(:,130:140,4),2);
snr>30
sum(snr>30)
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30&onOff>0,1:300,i)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30&onOff<0,1:300,i)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30&onOff>0,1:300,i)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30,1:300,i)')
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30,1:300,i)')
axis([0 300 -1 1])
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30&onOff<0,1:300,i)')
axis([0 300 -1 1])
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>30,1:300,i)')
axis([0 300 -1 1])
end
for i=1:9
subplot(3,3,i)
plot(meanFilter(snr>35,1:300,i)')
axis([0 300 -1 1])
end
filters(unit).name
unit
unit=1
filters(unit).name
load([filterspath,filters(unit).name])
LF=HighFilters';
plot(LF)
plot(LF')
plot(LF(37:45,:)')
plot(LF(28:36,:)')
HighSpikeCount(28:36)
clear
date='20130220_1';
unitspath=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/easy_formatted_units/'];
units=dir([unitspath,'*FFFlicker.mat']);
clear
date='20130220'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);
for unit=1:length(flashes)
length(flashes)
unit=1
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:82
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
trials=82;
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
kk=2;
for i=1:9
subplot(9,4,kk)
hold on
area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
line([500,500],[0,maxx],'Color','k','LineWidth',2);
line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
line([4500,4500],[0,maxx],'color','k','LineWidth',3);
line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
axis tight
set(gca,'XTick',0,'XTickLabel','')
title(['ND',int2str(i)])
kk=kk+4;
end
white_flash(:,1:4500)'
black_flash(:,1:4500)'
[white_flash(:,1:4500)' black_flash(:,1:4500)']
[white_flash(:,1:4500)'; black_flash(:,1:4500)']
res=zeros(46,9000,9);
for unit=1:length(flashes)
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
end
plot(res(:,:,4))
plot(res(:,:,4)')
plot(res(:,:,6)')
a=res(:,:,6)';
mean(a(1:400,:)
mean(a(1:400,:))
repmat(mean(a(1:400,:)),9000,1);
ans(1:100,)
ans(1:100,:)
ans(1:100,1)
a=a-b;
a=res(:,:,6)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
for i=1:9
subplot(3,3,i)
a=res(:,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
end
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(46+28,500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters';
cnt=1;
for i=1:9:80
meanFilter(unit+46,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
for unit=1:size(meanFilter,1)
if mean(meanFilter(unit,130:140,4),2)>0
a=reshape(max(meanFilter(unit,:,:)),9,1);
else
a=reshape(min(meanFilter(unit,:,:)),9,1);
end
tmp=reshape(meanFilter(unit,:,:),500,9)'./repmat(abs(a),1,500);
meanFilter(unit,:,:)=tmp';
end
snr=std(meanFilter(:,60:300,4),0,2)./std(meanFilter(:,1:20,4),0,2);
snr>30
date='20130220'
c=snr>30;
c=c(1:46);
i=1
subplot(3,3,i)
a=res(:,:,i)';
a=res(c,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
for i=1:9
subplot(3,3,i)
a=res(c,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
end
onOff=mean(meanFilter(:,130:140,4),2);
onOff
onOff=onOff(1:46);
for i=1:9
subplot(3,3,i)
a=res(c&onOff>1,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
end
for i=1:9
subplot(3,3,i)
a=res(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
end
figure
for i=1:9
subplot(3,3,i)
a=res(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
end
for i=1:9
subplot(3,3,i)
a=res(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i)])
end
for i=1:9
subplot(3,3,i)
a=res(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i)])
end
clear
date='20130220'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);
trials=82;
res=zeros(46+28,9000,9);
for unit=1:length(flashes)
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
end
date='20130220_1'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);
trials=81;
for unit=1:length(flashes)
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit+46,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
end
unit
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
for unit=1:length(flashes)
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=1:9:81
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit+46,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
end
clear
date='20130220'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(46+28,500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters';
cnt=1;
for i=1:9:80
meanFilter(unit+46,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
for unit=1:size(meanFilter,1)
if mean(meanFilter(unit,130:140,4),2)>0
a=reshape(max(meanFilter(unit,:,:)),9,1);
else
a=reshape(min(meanFilter(unit,:,:)),9,1);
end
tmp=reshape(meanFilter(unit,:,:),500,9)'./repmat(abs(a),1,500);
meanFilter(unit,:,:)=tmp';
end
snr=std(meanFilter(:,60:300,4),0,2)./std(meanFilter(:,1:20,4),0,2);
onOff=mean(meanFilter(:,130:140,4),2);
c=snr>30;
for i=1:9
subplot(3,3,i)
a=res(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i)])
end
res
clear
date='20130220'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);
trials=82;
res=zeros(46+28,9000,9);
for unit=1:length(flashes)
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
end
date='20130220_1'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);
trials=81;
for unit=1:length(flashes)
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=1:9:81
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit+46,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
end
date='20130220'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(46+28,500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters';
cnt=1;
for i=1:9:80
meanFilter(unit+46,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
for unit=1:size(meanFilter,1)
if mean(meanFilter(unit,130:140,4),2)>0
a=reshape(max(meanFilter(unit,:,:)),9,1);
else
a=reshape(min(meanFilter(unit,:,:)),9,1);
end
tmp=reshape(meanFilter(unit,:,:),500,9)'./repmat(abs(a),1,500);
meanFilter(unit,:,:)=tmp';
end
snr=std(meanFilter(:,60:300,4),0,2)./std(meanFilter(:,1:20,4),0,2);
onOff=mean(meanFilter(:,130:140,4),2);
c=snr>30;
for i=1:9
subplot(3,3,i)
a=res(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i)])
end
figure
for i=1:9
subplot(3,3,i)
a=res(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i)])
end
for i=1:9
subplot(3,3,i)
a=res(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i)])
end
get(gcf,'position')
figure
set(gcf,'position',[800    50   794   926])
get(gcf,'position')
figure
set(gcf,'position',[7    40   786   935])
for i=1:6
subplot(3,2,i)
a=res(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i)])
end
figure
set(gcf,'position',[7    40   786   935])
for i=1:6
subplot(3,2,i)
a=res(c&onOff>0,:,i+2)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2)])
end
close all
figure
set(gcf,'position',[800    50   794   926])
for i=1:6
subplot(3,2,i)
a=res(c&onOff<0,:,i+2)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2)])
end
figure
set(gcf,'position',[7    40   786   935])
for i=1:6
subplot(3,2,i)
a=res(c&onOff>0,:,i+2)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2)])
end
close all
figure
set(gcf,'position',[800    50   794   926])
for i=1:6
subplot(3,2,i)
a=res(c&onOff<0,:,i+2)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2),' n=',int2str(sum(c&onOff<0))])
end
figure
set(gcf,'position',[7    40   786   935])
for i=1:6
subplot(3,2,i)
a=res(c&onOff>0,:,i+2)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2),' n=',int2str(sum(c&onOff>0))])
end
clear
date='20130220'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*Moment*']);
trials=82;
res=zeros(46+28,5000,9);
for unit=1:length(flashes)
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit,:,:)=[white_flash(:,1:2500)'; black_flash(:,1:2500)'];
end
clear
date='20130220'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
moments=dir([path2data,'*Moment*']);
trials=82;
res=zeros(46+28,5000,9);
for unit=1:length(moments)
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit,:,:)=[white_flash(:,1:2500)'; black_flash(:,1:2500)'];
end
date='20130220_1'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
moments=dir([path2data,'*Moment*']);
trials=81;
for unit=1:length(moments)
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
res(unit+46,:,:)=[white_flash(:,1:2500)'; black_flash(:,1:2500)'];
end
date='20130220_1'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
moments=dir([path2data,'*Moment*']);
trials=81;
for unit=1:length(moments)
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=1:9:80
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit+46,:,:)=[white_flash(:,1:2500)'; black_flash(:,1:2500)'];
end
date='20130220'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(46+28,500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters';
cnt=1;
for i=1:9:80
meanFilter(unit+46,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
for unit=1:size(meanFilter,1)
if mean(meanFilter(unit,130:140,4),2)>0
a=reshape(max(meanFilter(unit,:,:)),9,1);
else
a=reshape(min(meanFilter(unit,:,:)),9,1);
end
tmp=reshape(meanFilter(unit,:,:),500,9)'./repmat(abs(a),1,500);
meanFilter(unit,:,:)=tmp';
end
snr=std(meanFilter(:,60:300,4),0,2)./std(meanFilter(:,1:20,4),0,2);
onOff=mean(meanFilter(:,130:140,4),2);
figure
set(gcf,'position',[800    50   794   926])
for i=1:6
subplot(3,2,i)
a=res(c&onOff<0,:,i+2)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2),' n=',int2str(sum(c&onOff<0))])
end
figure
set(gcf,'position',[7    40   786   935])
for i=1:6
subplot(3,2,i)
a=res(c&onOff>0,:,i+2)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2),' n=',int2str(sum(c&onOff>0))])
end
c=snr>30;
figure
set(gcf,'position',[800    50   794   926])
for i=1:6
subplot(3,2,i)
a=res(c&onOff<0,:,i+2)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2),' n=',int2str(sum(c&onOff<0))])
end
figure
set(gcf,'position',[7    40   786   935])
for i=1:6
subplot(3,2,i)
a=res(c&onOff>0,:,i+2)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2),' n=',int2str(sum(c&onOff>0))])
end
figure
set(gcf,'position',[800    50   794   926])
for i=1:6
subplot(3,2,i)
a=res(c&onOff<0,:,i+2)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2),' n=',int2str(sum(c&onOff<0))])
end
figure
set(gcf,'position',[7    40   786   935])
for i=1:6
subplot(3,2,i)
a=res(c&onOff>0,:,i+2)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i-2),' n=',int2str(sum(c&onOff>0))])
end
clear
date='20130220'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);
trials=82;
res=zeros(46+28,9000,9);
for unit=1:length(flashes)
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
end
date='20130220_1'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);
trials=81;
for unit=1:length(flashes)
% FULL FIELD FLASH
load([path2data,flashes(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,4501);
black_flash_tmp=zeros(trials,4501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
end
end
cnt=1;
for i=1:9:81
white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit+46,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
end
timing_flash=ones(1,9000);
diff(flips)
timing_flash=zeros(1,9000);
timing_flash(501:1984+501)=1;
timing_flash(5001:1984+5001)=-1;
plot(timing_flash)
flash
timing_flash=zeros(1,9000);
timing_flash(501:1984+501)=1;
timing_flash(5001:1984+5001)=-1;
flash=res;
clear res
date='20130220'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
moments=dir([path2data,'*Moment*']);
trials=82;
res=zeros(46+28,5000,9);
for unit=1:length(moments)
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=[1 11:9:81]
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit,:,:)=[white_flash(:,1:2500)'; black_flash(:,1:2500)'];
end
date='20130220_1'
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
moments=dir([path2data,'*Moment*']);
trials=81;
for unit=1:length(moments)
load([path2data,moments(unit).name])
clear white_flash black_flash
white_flash_tmp=zeros(trials,2501);
black_flash_tmp=zeros(trials,2501);
for trial=1:trials
spikes=cell2mat(spike_info.spike_times(trial));
flips=cell2mat(spike_info.flip_times(trial));
flips=flips(:,1);
conv=convolved(spikes,40,flips(end));
conv=conv(121:end-120);
for st=1:4:20
white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
end
end
cnt=1;
for i=1:9:80
white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
cnt=cnt+1;
end
white_flash=white_flash/5;
black_flash=black_flash/5;
res(unit+46,:,:)=[white_flash(:,1:2500)'; black_flash(:,1:2500)'];
end
moment
moments=res;
diff(flips)
timing_moments=zeros(1,5000);
timing_moments(501:17)=1;
timing_moments(2501:17+2501)=-1;
date='20130220'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(46+28,500,9);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters(:,[1:9 11:82])';
cnt=1;
for i=1:9:80
meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
for unit=1:length(filters)
load([filterspath,filters(unit).name])
LF=HighFilters';
cnt=1;
for i=1:9:80
meanFilter(unit+46,:,cnt)=mean(LF(i:i+8,:));
cnt=cnt+1;
end
end
for unit=1:size(meanFilter,1)
if mean(meanFilter(unit,130:140,4),2)>0
a=reshape(max(meanFilter(unit,:,:)),9,1);
else
a=reshape(min(meanFilter(unit,:,:)),9,1);
end
tmp=reshape(meanFilter(unit,:,:),500,9)'./repmat(abs(a),1,500);
meanFilter(unit,:,:)=tmp';
end
snr=std(meanFilter(:,60:300,4),0,2)./std(meanFilter(:,1:20,4),0,2);
onOff=mean(meanFilter(:,130:140,4),2);
save('S:\data\alexandra\MEA_data\flash_analysis\flash_moment.mat','flash','timing_flash','moments','timing_moments','meanFilter','snr','onOff')
clear
load('S:\data\alexandra\MEA_data\flash_analysis\flash_moment.mat']
clear
load('S:\data\alexandra\MEA_data\flash_analysis\flash_moment.mat')
c=snr>30;
figure
get(gcf,'position')
for i=3:8
subplot(3,2,i-2)
a=res(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
end
for i=3:8
subplot(3,2,i-2)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
end
close all
figure
set(gcf,'position',[ 1 31 1680 946])
for i=3:8
subplot(1,6,i-2)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
end
close all
figure
set(gcf,'position',[ 1 31 1680 946])
for i=3:8
subplot(6,1,i-2)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
end
close all
figure
set(gcf,'position',[ 1 31 1680 946])
for i=3:8
subplot(6,4,i-2)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
end
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
max(a)
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
i=2;
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
find(a(3315,:)<-1.5)
figure
i=2;
for j=1:15
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
figure
plot(a)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=moments(c&onOff<0,:,3)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
figure
plot(a)
find(a(3315,:)<-1.5)
a=flash(c&onOff>0,:,3)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
figure
plot(a)
find(a(3315,:)<-1.5)
a=flash(c&onOff>0,:,3)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
figure
plot(a)
find(a(5157,:)<-1.5)
c=snr>30;
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
a(:,7)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
a(:,7)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
a=flash(c&onOff>0,:,6)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
figure
plot(a)
find(a(5300,:)<-1.5)
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
a(:,[7 14])=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
a(:,[7 14])=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
a=moments(c&onOff>0,:,6)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
figure
plot(a)
find(a(664,:)<-1.5)
a=moments(c&onOff>0,:,6)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
figure
plot(a)
find(a(664,:)<-1.5)
a=moments(c&onOff>0,:,6)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
figure
plot(a)
find(a(6,:)<-1.5)
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,[7 14])=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
%     a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
%     a(:,[7 14])=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
a=moments(c&onOff>0,:,5)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
figure
a=moments(c&onOff>0,:,5)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
figure
plot(a)
a=moments(c&onOff>0,:,5)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
figure
plot(a)
find(a(3151,:)<-1.5)
find(a(3151,:)<-1)
figure
plot(a)
find(a(3151,:)<-0.96)
clear
load('S:\data\alexandra\MEA_data\flash_analysis\flash_moment.mat')
c=snr>30;
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
a(:,14)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
a(:,14)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
a=[a(4501:end,:), a(1:4500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
a(:,14)=[];
a=[a(2501:end,:), a(1:2500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=3;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%
a=[a(4501:end,:), a(1:4500,:)];
a=[a(4501:end,:); a(1:4500,:)];
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
a=[a(4501:end,:); a(1:4500,:)];
c=snr>30;
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
a=[a(4501:end,:); a(1:4500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
a(:,14)=[];
a=[a(2501:end,:); a(1:2500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
i=3;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
figure
imagesc(corr(a))
x=a;
x=cov(x');
x=cov(x);
[V,~]=eig(x);
pc_vectors=V(:,end-2:end);
pc1=spikes'*pc_vectors(:,end);
pc2=spikes'*pc_vectors(:,end-1);
pc1=a*pc_vectors(:,end);
pc2=a*pc_vectors(:,end-1);
figure
plot(pc1,pc2,'.')
pc1=a'*pc_vectors(:,end);
pc2=a'*pc_vectors(:,end-1);
x=a;
x=cov(x);
[V,~]=eig(x);
pc_vectors=V(:,end-2:end);
pc1=a*pc_vectors(:,end);
pc2=a'*pc_vectors(:,end-1);
figure
plot(pc1,pc2,'.')
figure
plot(pc_vectors)
i=3;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
figure
imagesc(corr(a))
i=3;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
figure
plot(a)
i=3;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
figure
plot(a)
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
figure
plot(a)
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
figure
plot(a)
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
figure
plot(a)
a1=max(a(500:1000,:))
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(a1./a2)
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(a1./a2)
hold on
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
plot(a2./a1,'r')
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
figure
plot(a)
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(a)
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
figure
plot(a)
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(a1./a2)
hold on
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(2480:2980,:));
a2=max(a(4930:5430,:));
plot(a2./a1,'r')
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(sort(a1./a2))
hold on
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(2480:2980,:));
a2=max(a(4930:5430,:));
plot(sort(a2./a1),'r')
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(sort(a1./a2))
hold on
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(2480:2980,:));
a2=max(a(4930:5430,:));
plot(sort(a2./a1),'r')
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(sort(a1./a2))
hold on
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(2480:2980,:));
a2=max(a(4930:5430,:));
plot(sort(a1./a2),'r')
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
figure
plot(a)
figure
plot(sort(a1./a2))
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(sort(a1./a2))
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
figure
plot(a)
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(sort(a1./a2))
hold on
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(2480:2980,:));
a2=max(a(4930:5430,:));
plot(sort(a2./a2),'r')
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(sort(a1./a2))
hold on
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(2480:2980,:));
a2=max(a(4930:5430,:));
plot(sort(a2./a2),'r')
i=5;
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
figure
plot(sort(a1./a2))
hold on
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(2480:2980,:));
a2=max(a(4930:5430,:));
plot(sort(a2./a1),'r')
figure
for i=3:8
subplot(2,3,i-2)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
plot(sort(a1./a2))
hold on
i=5;
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(2480:2980,:));
a2=max(a(4930:5430,:));
plot(sort(a2./a1),'r')
title(['ND',int2str(10-i)])
end
figure
for i=3:8
subplot(2,3,i-2)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(500:1000,:));
a2=max(a(7000:7500,:));
plot(sort(a1./a2))
hold on
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
% b=repmat(max(a),9000,1);
% a=a./abs(b);
a1=max(a(2480:2980,:));
a2=max(a(4930:5430,:));
plot(sort(a2./a1),'r')
title(['ND',int2str(10-i)])
end
c=snr>25;
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
a=[a(4501:end,:); a(1:4500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
a(:,14)=[];
a=[a(2501:end,:); a(1:2500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
c=snr>35;
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
a=[a(4501:end,:); a(1:4500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
a(:,14)=[];
a=[a(2501:end,:); a(1:2500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
c=snr>35;
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
a=[a(4501:end,:); a(1:4500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
a=[a(2501:end,:); a(1:2500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
c=snr>30;
close all
figure
set(gcf,'position',[ 1 31 1680 946])
cnt=1;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,8)=[];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=3;
for i=3:8
subplot(6,4,cnt)
a=flash(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
a=[a(4501:end,:); a(1:4500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
cnt=2;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff<0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
cnt=cnt+4;
end
cnt=4;
for i=3:8
subplot(6,4,cnt)
a=moments(c&onOff>0,:,i)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
a=[a(2501:end,:); a(1:2500,:)];
plot(a)
axis tight
title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
cnt=cnt+4;
end
clear
flag=1;
last_copied_mcd=1;
last_copied_heka=1;
mcd_orig_path='\\Cin-tm-1ephys\Alexandra\20130224\MEA_mcd\';
heka_orig_path='\\Cin-tm-1ephys\Alexandra\20130224\HEKA\';
mcd_dest_path='S:\data\alexandra\MEA_data\20130224\MEA_mcd\';
heka_dest_path='S:\data\alexandra\MEA_data\20130224\HEKA\';
if ~exist(mcd_dest_path,'file')
mkdir(mcd_dest_path);
end
if ~exist(heka_dest_path,'file')
mkdir(heka_dest_path);
end
while flag==1
mcd=dir([mcd_orig_path,'*.mcd']);
heka=dir([heka_orig_path,'*.phys']);
a=length(mcd);
b=length(heka);
while last_copied_heka<b
copyfile([heka_orig_path,heka(last_copied_heka).name],[heka_dest_path,heka(last_copied_heka).name])
last_copied_heka=last_copied_heka+1;
end
while last_copied_mcd<a
copyfile([mcd_orig_path,mcd(last_copied_mcd).name],[mcd_dest_path,mcd(last_copied_mcd).name])
last_copied_mcd=last_copied_mcd+1;
end
k=clock;
fprintf('Next Step\tHEKA=%d\tMCD=%d\t\t%d:%d\n',b,a,k(4),k(5))
pause(60)
end
clc
clear