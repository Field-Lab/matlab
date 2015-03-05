%% Linear Filter, full field.

%get linear filter
clear
date='20120830';
unitspath=['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\'];
units=dir([unitspath,'*FFFlicker.mat']);

path2save=['S:\user\alexandra\MEA_data\',date,'\FFFlicker_LF\'];

if ~exist(path2save,'file')
    mkdir(path2save);
end
load([unitspath,units(1).name])
filter_length=500; %filter length in ms

for units_cnt=1:length(units)
    units(units_cnt).name
    if ~exist([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter.mat'],'file')
        load([unitspath,units(units_cnt).name])
        trials=size(spike_info.flip_times,1);
        
        HighSpikeCount=zeros(1,trials);
        HighFilters=zeros(filter_length,trials);
        cntH=1;
        
        for trial_counter=1:trials
            flips=spike_info.flip_times{trial_counter}(1:29*60,1); % flips in ms, rounded
            pxls=(((spike_info.flip_times{trial_counter}(1:29*60,2))-30)/30); % pxs in ms, normalized -1:1
            spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
            spikes(spikes<filter_length|spikes>flips(end))=[];
            if ~isempty(spikes)
                tmp=zeros(flips(end)+1,2); % 1 column - color, 2 column - spikes
                tmp([1; flips(1:end-1)],1)=pxls;
                
                % fill stimulus down to ms with brightness values where spikes
                % happen
                for i=1:50
                    subscr=(flips(1:end-1)+i);
                    a=flips(1:end-1)+i;
                    while a(end)>length(tmp)
                        a(end)=[];
                    end
                    subscr(tmp(a,1)~=0)=[];
                    tmp(subscr,1)=tmp(subscr-1,1);
                end
                
                tmp(spikes,2)=1;
                
                n=zeros(length(spikes),filter_length); % filter
                for i=1:filter_length
                    n(:,i)=tmp(spikes-i+1,1);
                    %     n(i)=sum(tmp(c-i+1,1))/all_spikes;
                end
                n=[spikes, n];
                tmp=n(n(:,1)<flips(1740)&n(:,1)>filter_length,2:end);
                HighFilters(1:filter_length,trial_counter)=sum(tmp);
                HighSpikeCount(trial_counter)=size(tmp,1);
            end
        end
        save([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter.mat'],'HighFilters','HighSpikeCount')
    end
end


% plotting
cd('S:\user\alexandra\scripts')
clear
date='20120329';
filter_length=500;
nd='876543212345';
filterspath=['F:\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);

path2save=['F:\',date,'\FFFlicker_LF_plots\'];

if ~exist(path2save,'file')
    mkdir(path2save);
end

close(figure(1))
figure(1)
set(gcf,'Position', [1 31 1680 946])
[rows,cols]=opt_subplots(length(nd));
% simple filter plot
for filters_cnt=1:length(filters)
    if ~exist([path2save,filters(filters_cnt).name(1:end-3),'emf'],'file')
        load([filterspath,filters(filters_cnt).name])
        
        minn=min(min(HighFilters(:)));
        maxx=max(max(HighFilters(:)));
        cnt=1;
        for i=1:4:length(nd)*4
            k=subplot(rows,cols,cnt);
            goto=min(i+3,size(HighFilters,2));
            plot(HighFilters(:,i:goto))
            line([0 500],[0 0],'color',[1 1 1]*0.5)
            hleg=legend(int2str(HighSpikeCount(i:goto)'));
            legend('boxoff')
            set(hleg,'Fontsize',8)
            title(['ND',nd(cnt)])
            %         axis tight
            %         axis([0 500 -100 150])
            axis([0 500 minn maxx])
            cnt=cnt+1;
        end
        
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title([filters(filters_cnt).name(1:end-4),' ',int2str(filter_length),'ms   High CONTRAST 30s'],'FontSize',15,'FontWeight','bold','Interpreter','none')
        
        saveas(gcf,[path2save,filters(filters_cnt).name(1:end-3),'emf'])
    end
end


