%% Linear Filter, full field.

%get linear filter
clear
date='20120908_1';
unitspath=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/easy_formatted_units/'];
units=dir([unitspath,'*FFFlicker.mat']);

path2save=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/FFFlicker_LF/'];

if ~exist(path2save,'file')
    mkdir(path2save);
end

filter_length=500; %filter length in ms

for units_cnt=1:length(units)
    
    load([unitspath,units(units_cnt).name])
    trials=size(spike_info.flip_times,1);
    
    HighSpikeCount=zeros(1,trials);
    LowSpikeCount=zeros(1,trials);
    
    for trial_counter=1:trials
        flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
        pxls=(((spike_info.flip_times{trial_counter}(:,2))-30)/30); % pxs in ms, normalized -1:1
        spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
        
        spikes(spikes<filter_length|spikes>flips(end))=[];
        
        tmp=zeros(flips(end)+1,2); % 1 column - color, 2 column - spikes
        tmp([1; flips(1:end-1)],1)=pxls;
        
        % fill stimulus down to ms with brightness values where spikes
        % happen
        for i=1:50
            subscr=(flips(1:end-1)+i);
            subscr(tmp(flips(1:end-1)+i,1)~=0)=[];
            tmp(subscr)=tmp(subscr-1);
        end
        
        tmp(spikes,2)=1;
        
        n=zeros(length(spikes),filter_length); % filter
        for i=1:filter_length
            n(:,i)=tmp(spikes-i+1,1);
            %     n(i)=sum(tmp(c-i+1,1))/all_spikes;
        end
        n=[spikes, n];
        % plot(b)
        contrast_changes=flips(1:1800:end);
        tmp=0;
        HighFilter=0;
        if length(flips)>17000
            st=10;
        else
            st=8;
        end
        for i=1:2:st
            tmp=n(n(:,1)<contrast_changes(i+1)&n(:,1)>contrast_changes(i)+filter_length,2:end);
            HighFilter=HighFilter+sum(tmp);
            HighSpikeCount(trial_counter)=HighSpikeCount(trial_counter)+size(tmp,1);
        end
        
        LowFilter=0;
        for i=2:2:st
            tmp=n(n(:,1)<contrast_changes(i+1)&n(:,1)>contrast_changes(i)+filter_length,2:end);
            LowFilter=LowFilter+sum(tmp);
            LowSpikeCount(trial_counter)=LowSpikeCount(trial_counter)+size(tmp,1);
        end
        
        %         HighFilter=HighFilter(end:-1:1)/abs(min(HighFilter)); % mean kernel
        HighFilters(1:filter_length,trial_counter)=HighFilter;
        %         LowFilter=LowFilter(end:-1:1)/abs(min(LowFilter)); % mean kernel
        LowFilters(1:filter_length,trial_counter)=LowFilter;
    end
    save([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter.mat'],'HighFilters','HighSpikeCount','LowFilters','LowSpikeCount')
    
end



% plotting
clear
date='20120908_1';
filter_length=500;
nd='43212345';
filterspath=['F:\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);

path2save=['F:\',date,'\FFFlicker_LF_plots\'];

if ~exist(path2save,'file')
    mkdir(path2save);
end

% simple filter plot
for filters_cnt=1:length(filters)
    load([filterspath,filters(filters_cnt).name])
    trials=size(HighFilters,2);
    close(figure(1))
    figure(1)
    set(gcf,'Position', [1 31 1680 946])
    minn=min(HighFilters(:));
    maxx=max(HighFilters(:));
    cnt=1;
    for i=1:4:trials
        subplot(2,4,cnt)
        goto=min(i+3,trials);
        plot(HighFilters(:,i:goto))
        legend(int2str(HighSpikeCount(i:goto)'))
        title(['ND',nd(cnt)])
        axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([filters(filters_cnt).name(1:end-4),' ',int2str(filter_length),'ms'],'FontSize',15,'FontWeight','bold','Interpreter','none')
    
    saveas(gcf,[path2save,filters(filters_cnt).name(1:end-3),'emf'])
    
end






% special for 08092012   simple filter plot
path2save=['F:\',date,'\FFFlicker_LF_plots\'];
filter_length=300;
for filters_cnt=1:length(filters)
    load([filterspath,filters(filters_cnt).name])
    trials=size(HighFilters,2);
    close(figure(1))
    figure(1)
    set(gcf,'Position', [1 31 1680 946])
    minn=min(HighFilters(:));
    maxx=max(HighFilters(:));
    cnt=1;
    for i=1:4:12
        subplot(2,5,cnt)
        plot(HighFilters(:,i:i+3))
        legend(int2str(HighSpikeCount(i:i+3)'))
        title(['ND4 cycle ',int2str(cnt)])
        axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end

    subplot(2,5,4)
    plot(HighFilters(:,13))
    legend(int2str(HighSpikeCount(13)'))
    title(['ND3 single'])
    axis([0 filter_length minn maxx])

    cnt=5;
    for i=14:4:22
        subplot(2,5,cnt)
        plot(HighFilters(:,i:i+3))
        legend(int2str(HighSpikeCount(i:i+3)'))
        title(['ND',nd(cnt-4)])
        axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end

    subplot(2,5,8)
    plot(HighFilters(:,26:27))
    legend(int2str(HighSpikeCount(26:27)'))
    title(['ND2, 3 empty'])
    axis([0 filter_length minn maxx])
    
   
    cnt=9;
    for i=28:4:32
        subplot(2,5,cnt)
        plot(HighFilters(:,i:i+3))
        legend(int2str(HighSpikeCount(i:i+3)'))
        title('ND3')
        axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end
    
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([filters(filters_cnt).name(1:end-4),' ',int2str(filter_length),'ms'],'FontSize',15,'FontWeight','bold','Interpreter','none')
    
    saveas(gcf,[path2save,filters(filters_cnt).name(1:end-3),'emf'])
    
end


path2save=['F:\',date,'\FFFlicker_LF_plots_trial_group\'];
filter_length=300;
for filters_cnt=1:length(filters)
    load([filterspath,filters(filters_cnt).name])
    trials=size(HighFilters,2);
    close(figure(1))
    figure(1)
    set(gcf,'Position', [1 31 1680 946])
    minn=min(HighFilters(:));
    maxx=max(HighFilters(:));
    cnt=1;
    for i=1:4
        subplot(2,4,cnt)
        plot(HighFilters(:,[i:4:i+8,i+13]))
        hold on
        plot(HighFilters(:,i+17),'k','linewidth',2)
        plot(HighFilters(:,i+21),'m','linewidth',2)
        legend(int2str(HighSpikeCount([i:4:i+8,i+13,i+17,i+21])'))
        title(['ND4 4 times, ND3, ND2, trial ',int2str(cnt),' '])
        axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end

    subplot(2,4,5)
    plot(HighFilters(:,13))
    legend(int2str(HighSpikeCount(13)'))
    title(['ND3 single'])
    axis([0 filter_length minn maxx])

    subplot(2,4,6)
    plot(HighFilters(:,26:27))
    legend(int2str(HighSpikeCount(26:27)'))
    title(['ND2, 3 empty'])
    axis([0 filter_length minn maxx])
    
   
    cnt=7;
    for i=28:4:32
        subplot(2,4,cnt)
        plot(HighFilters(:,i:i+3))
        legend(int2str(HighSpikeCount(i:i+3)'))
        title('ND3')
        axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end
    
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([filters(filters_cnt).name(1:end-4),' ',int2str(filter_length),'ms'],'FontSize',15,'FontWeight','bold','Interpreter','none')
    
    saveas(gcf,[path2save,filters(filters_cnt).name(1:end-3),'emf'])
    
end








%% LF from partial sequences

%get linear filter
clear
date='20120908_1';
unitspath=['S:\user\alexandra\MEA_data\',date,'/easy_formatted_units/'];
units=dir([unitspath,'*FFFlicker.mat']);

path2save=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/FFFlicker_LF/'];

if ~exist(path2save,'file')
    mkdir(path2save);
end

filter_length=500; %filter length in ms

units_cnt=3;
load([unitspath,units(units_cnt).name])
trials=size(spike_info.flip_times,1);
trials=2;
HighSpikeCount=zeros(1,trials*4);
LowSpikeCount=zeros(1,trials*4);
cnt=1;cnt1=1;
for trial_counter=26:27
    flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
    pxls=(((spike_info.flip_times{trial_counter}(:,2))-30)/30); % pxs in ms, normalized -1:1
    spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
    
    spikes(spikes<filter_length|spikes>flips(end))=[];
    
    tmp=zeros(flips(end)+1,2); % 1 column - color, 2 column - spikes
    tmp([1; flips(1:end-1)],1)=pxls;
    
    % fill stimulus down to ms with brightness values where spikes
    % happen
    for i=1:50
        subscr=(flips(1:end-1)+i);
        subscr(tmp(flips(1:end-1)+i,1)~=0)=[];
        tmp(subscr)=tmp(subscr-1);
    end
    
    tmp(spikes,2)=1;
    
    n=zeros(length(spikes),filter_length); % filter
    for i=1:filter_length
        n(:,i)=tmp(spikes-i+1,1);
        %     n(i)=sum(tmp(c-i+1,1))/all_spikes;
    end
    n=[spikes, n];
    % plot(b)
    contrast_changes=flips(1:1800:end);
    tmp=0;
    for i=1:2:8
        tmp=n(n(:,1)<contrast_changes(i+1)&n(:,1)>contrast_changes(i)+filter_length,2:end);
        HighFilters(1:filter_length,cnt)=sum(tmp);
        HighSpikeCount(cnt)=size(tmp,1);
        cnt=cnt+1;
    end
    for i=2:2:8
        tmp=n(n(:,1)<contrast_changes(i+1)&n(:,1)>contrast_changes(i)+filter_length,2:end);
        LowFilters(1:filter_length,cnt1)=sum(tmp);
        LowSpikeCount(cnt1)=size(tmp,1);
        cnt1=cnt1+1;
    end
    
end

subplot(2,2,1)
plot(HighFilters(:,1:4),'linewidth',2)
legend(int2str(HighSpikeCount(1:4)'))
title([units(units_cnt).name(1:27),', High contrast 1 trial, ND2, seq 2 1 empty 4'],'interpreter','none')
subplot(2,2,2)
plot(LowFilters(:,1:4),'linewidth',2)
legend(int2str(LowSpikeCount(1:4)'))
title([units(units_cnt).name(1:27),', Low contrast 1 trial, ND2, seq 2 1 empty 4'],'interpreter','none')
subplot(2,2,3)
plot(HighFilters(:,5:8),'linewidth',2)
legend(int2str(HighSpikeCount(5:8)'))
title([units(units_cnt).name(1:27),', High contrast 2 trial, ND2, seq 2 1 empty 4'],'interpreter','none')
subplot(2,2,4)
plot(LowFilters(:,5:8),'linewidth',2)
legend(int2str(LowSpikeCount(5:8)'))
title([units(units_cnt).name(1:27),', Low contrast 2 trial, ND2, seq 2 1 empty 4'],'interpreter','none')


% simple filter plot
% 7 ND, 6 trials in each, 4 dif. sequences in each
for nd=[3,6]
    
    figure
    set(gcf,'Position', [1 31 1680 946])
    minn=min(HighFilters(:));
    maxx=max(HighFilters(:));
    cnt=1;
    tr='123456';
    
    st=8-nd;
    for i=24*(st-1)+1:4:24*st
        subplot(2,4,cnt)
        plot(HighFilters(:,i:i+3))
        legend(int2str(HighSpikeCount(i:i+3)'))
        hold on
        plot(LowFilters(:,i:i+3),'linewidth',2)
        title(['ND',int2str(nd),', trial ',tr(cnt)])
        axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end
    subplot(2,4,7)
    a=[];
    for i=24*(st-1)+1:4:24*st
        a=[a sum(HighFilters(:,i:i+3),2)];
    end
    plot(a)
    title(['ND',int2str(nd),', 6 trials High contrast'])
    subplot(2,4,8)
    a=[];
    for i=24*(st-1)+1:4:24*st
        a=[a sum(LowFilters(:,i:i+3),2)];
    end
    plot(a)
    title(['ND',int2str(nd),', 6 trials Low contrast'])
    
    
    figure
    set(gcf,'Position', [1 31 1680 946])
    minn=min(HighFilters(:));
    maxx=max(HighFilters(:));
    cnt=1;
    seq='1234';
    st=8-nd;
    for i=24*(st-1)+1:24*(st-1)+4
        subplot(2,3,cnt)
        plot(HighFilters(:,i:4:24*st))
        legend(int2str(HighSpikeCount(i:4:24*st)'))
        hold on
        plot(LowFilters(:,i:4:24*st),'linewidth',2)
        title(['ND',int2str(nd),', seq ',seq(cnt)])
        axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end
    
    subplot(2,3,5)
    a=[];
    for i=24*(st-1)+1:24*(st-1)+4
        a=[a sum(HighFilters(:,i:4:24*st),2)];
    end
    plot(a)
    title(['ND',int2str(nd),', 4 sequences High contrast'])
    
    subplot(2,3,6)
    a=[];
    for i=24*(st-1)+1:24*(st-1)+4
        a=[a sum(LowFilters(:,i:4:24*st),2)];
    end
    plot(a)
    title(['ND',int2str(nd),', 4 sequences Low contrast'])
end




% simple filter plot
% 7 ND, 6 trials in each, 4 dif. sequences in each, taken 5-15s and 15-25s
% intervals from each sequence
for nd=[3,6]
    
    figure
    set(gcf,'Position', [1 31 1680 946])
    minn=min(HighFilters(:));
    maxx=max(HighFilters(:));
    cnt=1;
    tr='123456';
    
    st=8-nd;
    for i=48*(st-1)+1:8:48*st
        subplot(2,4,cnt)
        plot(HighFilters(:,i:i+3))
        legend(int2str(HighSpikeCount(i:i+3)'))
%         hold on
%         plot(LowFilters(:,i:i+3),'linewidth',2)
%         title(['ND',int2str(nd),', trial ',tr(cnt)])
%         axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end
    subplot(2,4,7)
    a=[];
    for i=24*(st-1)+1:4:24*st
        a=[a sum(HighFilters(:,i:i+3),2)];
    end
    plot(a)
    title(['ND',int2str(nd),', 6 trials High contrast'])
    subplot(2,4,8)
    a=[];
    for i=24*(st-1)+1:4:24*st
        a=[a sum(LowFilters(:,i:i+3),2)];
    end
    plot(a)
    title(['ND',int2str(nd),', 6 trials Low contrast'])
    
    
    figure
    set(gcf,'Position', [1 31 1680 946])
    minn=min(HighFilters(:));
    maxx=max(HighFilters(:));
    cnt=1;
    seq='1234';
    st=8-nd;
    for i=24*(st-1)+1:24*(st-1)+4
        subplot(2,3,cnt)
        plot(HighFilters(:,i:4:24*st))
        legend(int2str(HighSpikeCount(i:4:24*st)'))
        hold on
        plot(LowFilters(:,i:4:24*st),'linewidth',2)
        title(['ND',int2str(nd),', seq ',seq(cnt)])
        axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end
    
    subplot(2,3,5)
    a=[];
    for i=24*(st-1)+1:24*(st-1)+4
        a=[a sum(HighFilters(:,i:4:24*st),2)];
    end
    plot(a)
    title(['ND',int2str(nd),', 4 sequences High contrast'])
    
    subplot(2,3,6)
    a=[];
    for i=24*(st-1)+1:24*(st-1)+4
        a=[a sum(LowFilters(:,i:4:24*st),2)];
    end
    plot(a)
    title(['ND',int2str(nd),', 4 sequences Low contrast'])
end










%% High-Low contrast
% plotting
clear
date='20120714';
filter_length=500;
nd='76543212345';
filterspath=['F:\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);

path2save=['F:\',date,'\FFFlicker_LF_plots\'];

if ~exist(path2save,'file')
    mkdir(path2save);
end

% simple filter plot
for filters_cnt=1:length(filters)
    load([filterspath,filters(filters_cnt).name])
    trials=size(HighFilters,2);
    close(figure(1))
    figure(1)
    set(gcf,'Position', [1 31 1680 946])
    minn=min(HighFilters(:));
    maxx=max(HighFilters(:));
    cnt=1;
    for i=1:6:trials
        subplot(2,4,cnt)
        pl=mean(HighFilters(:,i:i+5),2);
        pl=pl/max(pl);
        plot(pl)
        hold on
        pl=mean(LowFilters(:,i:i+5),2);
        pl=pl/max(pl);
        plot(pl,'r')
%         legend(int2str(HighSpikeCount(i:goto)'),int2str(LowSpikeCount(i:goto)'))
        title(['ND',nd(cnt)])
%         axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([filters(filters_cnt).name(1:end-4),' ',int2str(filter_length),'ms'],'FontSize',15,'FontWeight','bold','Interpreter','none')
    
    saveas(gcf,[path2save,filters(filters_cnt).name(1:end-3),'emf'])
    
end
