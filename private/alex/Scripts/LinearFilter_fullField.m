%% COMPARE DIFFERENT ND LINEAR PARTS

%% Linear part

%getting linear filter
clear
date='20120411';
unitspath=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/easy_formatted_units/'];
units=dir([unitspath,'*FFFlicker.mat']);

filter_length=500;
trialsToTake=49;

for units_cnt=1:length(units)
    
    load([unitspath,units(units_cnt).name])
    
    HighSpikeCount=zeros(1,trialsToTake);
    LowSpikeCount=zeros(1,trialsToTake);
    
    for trial_counter=1:trialsToTake
        flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
        pxls=(((spike_info.flip_times{trial_counter}(:,2))-30)/30); % pxs in ms, normalized -1:1
        % b=((spike_info.flip_times{trial_counter}(:,2))); % pxs in ms, not normalized
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
        for i=1:2:8
            tmp=n(n(:,1)<contrast_changes(i+1)&n(:,1)>contrast_changes(i)+filter_length,2:end);
            HighFilter=HighFilter+sum(tmp);
            HighSpikeCount(trial_counter)=HighSpikeCount(trial_counter)+size(tmp,1);
        end
        
        LowFilter=0;
        for i=2:2:8
            tmp=n(n(:,1)<contrast_changes(i+1)&n(:,1)>contrast_changes(i)+filter_length,2:end);
            LowFilter=LowFilter+sum(tmp);
            LowSpikeCount(trial_counter)=LowSpikeCount(trial_counter)+size(tmp,1);
        end
        
        HighFilters(1:filter_length,trial_counter)=HighFilter;
        LowFilters(1:filter_length,trial_counter)=LowFilter;
    end
    save(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/FFFlicker_LF/',units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter.mat'],'HighFilters','HighSpikeCount','LowFilters','LowSpikeCount')
end

% plotting
clear
date='20120627';
filter_length=500;

filterspath=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/FFFlicker_LF/'];
filters=dir([filterspath,'*filter.mat']);

nd='87654321';
% simple filter plot
for filters_cnt=1:length(filters)
    load([filterspath,filters(filters_cnt).name])
    close(figure(1))
    figure(1)
    set(gcf, 'Units','normalized','Position',[0.1    0.1    0.85    0.75])
    minn=min(HighFilters(:));
    maxx=max(HighFilters(:));
    cnt=1;
    for i=1:6:48
        subplot(2,4,cnt)
        plot(HighFilters(:,i:i+5))
        legend(int2str(HighSpikeCount(i:i+5)'))
        title(['ND',nd(cnt)])
        axis([0 filter_length minn maxx])
        cnt=cnt+1;
    end
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([filters(filters_cnt).name(1:end-4),' ',int2str(filter_length),'ms'])    
    saveas(gcf,['F:\',date,'\processing\filter_jitter_ND_pics_',int2str(filter_length),'\',filters(filters_cnt).name(1:end-3),'emf'])
    
end



% normalized within nd by the peak (compare timing)
for filters_cnt=1:length(filters)
    load([filterspath,filters(filters_cnt).name])
    close(figure(1))
    figure(1)
    set(gcf,'Position', [1 31 1680 946])    
    cnt=1;
    for i=1:4:36
        subplot(3,3,cnt)
        a=min(min(HighFilters(:,i:i+3)));        
        for j=i:i+3
            HighFilters(HighFilters(1:300,j)<0,j)=HighFilters(HighFilters(1:300,j)<0,j)*(a/min(HighFilters(:,j)));            
        end
        a=max(max(HighFilters(:,i:i+3)));        
        for j=i:i+3
            HighFilters(HighFilters(1:300,j)>0,j)=HighFilters(HighFilters(1:300,j)>0,j)*(a/max(HighFilters(:,j)));            
        end
        plot(HighFilters(:,i:i+3),'LineWidth',2.5)
        %     legend(['1234']')
        legend(int2str(HighSpikeCount(i:i+3)'))
        if i~=33
            title(['ND',int2str(9-cnt)])
        else
            title('ND2')
        end
        axis tight
        cnt=cnt+1;
    end
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([filters(filters_cnt).name(1:end-4),' ',int2str(filter_length),'ms, adjusted min and max (0-300ms) to compare timing'],'FontSize',15,'FontWeight','bold','Interpreter','none')
    
    saveas(gcf,['F:\',date,'\processing\filter_jitter_ND_pics_',int2str(filter_length),'_timing\',filters(filters_cnt).name(1:end-3),'emf'])
    
end



% fitting
% parameters: a1 - amplitude, b1 - delay (position of the peak), c1 - width
clear
date='20120301';
filter_length=500;
filterspath=['F:\',date,'\processing\linear_filters_raw_normal_500\'];
filters=dir([filterspath,'*filter.mat']);
on=0;off=0;
intervals=100;
trials=40;
lower_border=-10000;
upper_border=10000;
spikeLowerBorder=50;
peaks_inds=zeros(length(filters),trials);
a1=zeros(length(filters),trials);
b1=zeros(length(filters),trials);
c1=zeros(length(filters),trials);
on_off=zeros(length(filters),1)-10;
names=cell(1,length(filters));
for filters_cnt=1:length(filters)
    filters_cnt
    load([filterspath,filters(filters_cnt).name])
    names{filters_cnt}=filters(filters_cnt).name;
    [~, ind_max]=max(HighFilters(1:200,:));
    [~, ind_min]=min(HighFilters(1:200,:));
    if sum(10*std(HighFilters(1:50,5:20))<std(HighFilters(100:200,5:20)))
        if sum(ind_max>ind_min)>25
            polarity='OFF';% off
            on_off(filters_cnt)=-1;
            off=off+1;
            [~,ind_peak]=min(HighFilters(50:200,1:trials));
            peaks_inds(filters_cnt,1:trials)=ind_peak+49;
            for trial=1:trials
                if HighSpikeCount(trial)>spikeLowerBorder
                    x=(max(peaks_inds(filters_cnt,trial)-intervals,1):peaks_inds(filters_cnt,trial)+intervals)';
                    for_fit=-HighFilters(x,trial);
                    fit_res=fit(x,for_fit,'gauss1','Lower',[lower_border,lower_border,lower_border],'Upper',[upper_border,upper_border,upper_border]);
                    a1(filters_cnt,trial)=-fit_res.a1;
                    b1(filters_cnt,trial)=fit_res.b1;
                    c1(filters_cnt,trial)=fit_res.c1;
                end
            end
        elseif sum(ind_max>ind_min)<15
            polarity='ON';% on
            on_off(filters_cnt)=1;
            on=on+1;
            [~,ind_peak]=max(HighFilters(50:200,1:trials));
            peaks_inds(filters_cnt,1:trials)=ind_peak+49;
            for trial=1:trials
                if HighSpikeCount(trial)>spikeLowerBorder
                    x=(max(peaks_inds(filters_cnt,trial)-intervals,1):peaks_inds(filters_cnt,trial)+intervals)';
                    for_fit=HighFilters(x,trial);
                    fit_res=fit(x,for_fit,'gauss1','Lower',[lower_border,lower_border,lower_border],'Upper',[upper_border,upper_border,upper_border]);
                    a1(filters_cnt,trial)=fit_res.a1;
                    b1(filters_cnt,trial)=fit_res.b1;
                    c1(filters_cnt,trial)=fit_res.c1;
                end
            end
        else
            polarity='??';% not clear
            on_off(filters_cnt)=0;
        end
    end
    
end

save(['F:\',date,'\processing\linear_filters_fit_v1\gaussian_parameters.mat'],'a1','b1','c1','names','peaks_inds')

% plotting fit
load(['F:\',date,'\processing\linear_filters_fit_v1\gaussian_parameters.mat'])
for filters_cnt=1:length(filters)
    if abs(on_off(filters_cnt))==1
        load([filterspath,filters(filters_cnt).name])
        close (figure(1))
        figure(1)
        set(gcf,'Position', [1 31 1680 946])
        for trial=1:40
            subplot(7,6,trial)
            if HighSpikeCount(trial)>spikeLowerBorder
                col='b';
            else
                col='k';
            end
            plot(HighFilters(:,trial),col,'LineWidth',2)
            if a1(filters_cnt,trial)~=0 && b1(filters_cnt,trial)~=0 && c1(filters_cnt,trial)~=0
                hold on
                x=(peaks_inds(filters_cnt,trial)-intervals:peaks_inds(filters_cnt,trial)+intervals)';
                y=a1(filters_cnt,trial)*exp(-((x-b1(filters_cnt,trial))/c1(filters_cnt,trial)).^2);
                plot(x,y,'r','LineWidth',2)
            end
            axis tight
            title(['ND', int2str(9-floor((trial+3)/4)),'  trial ',int2str(mod(trial+3,4)+1)])
        end
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title([filters(filters_cnt).name, '    GAUSSIAN FIT of High contrast filter; blue: data, red: fit, black: too few spikes'],'Interpreter','None','FontSize',15,'FontWeight','bold')
        saveas(gcf,['F:\',date,'\processing\linear_filters_fit_v1\',filters(filters_cnt).name(1:end-3),'emf'])
    end
end


% plot delays
date='20120330';
load(['F:\',date,'\processing\linear_filters_fit\gaussian_parameters.mat'])
for i=1:size(b1,1)
    for nd=1:4:size(b1,2)-4
        plot(nd:nd+3,b1(i,nd:nd+3),'-*')
        hold on
    end
end


for i=1:size(b1,1)
    for nd=1:4:size(b1,2)-4
        if nnz(b1(i,nd:nd+3))==4
            plot(nd:nd+3,b1(i,nd:nd+3),'-*')
            hold on
        end
    end
end


% crosscorrelations

clear
date='20120329';
unitspath=['F:\',date,'\processing\easy_formatted_units\'];
units=dir([unitspath,'*FFFlicker.mat']);

% binary spike events
for units_cnt=1:length(units)
    
    load([unitspath,units(units_cnt).name])
       
    cnt=1;
    close (figure(1))
    figure(1)
    set(gcf,'Position', [1 31 1680 946])
    for i=1:4:32
        flips=spike_info.flip_times{i}(1:1800:end,1); % flips in ms, rounded
        trial1=i;
        trial2=i+1;
        tmp=spike_info.spike_times{trial1};
        spikes1=zeros(1,30000);
        for j=2:2:length(flips)
            spikes1(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999)=spikes1(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999)+1;
        end
        tmp=spike_info.spike_times{trial2};
        spikes2=zeros(1,30000);
        for j=2:2:length(flips)
            spikes2(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999)=spikes2(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999)+1;
        end

        [c_ww,lags] = xcorr(spikes1,spikes2,250,'coef');
        subplot(2,4,cnt)
        plot(lags,c_ww)
        hold on
        trial1=i;
        trial2=i+3;
        tmp=spike_info.spike_times{trial1};
        spikes1=zeros(1,30000);
        for j=2:2:length(flips)
            spikes1(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999)=spikes1(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999)+1;
        end
        tmp=spike_info.spike_times{trial2};
        spikes2=zeros(1,30000);
        for j=2:2:length(flips)
            spikes2(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999)=spikes2(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999)+1;
        end

        [c_ww,lags] = xcorr(spikes1,spikes2,250,'coef');
        subplot(2,4,cnt)
        plot(lags,c_ww,'r')
        title(['ND', int2str(9-cnt)])
        if cnt==1
            legend(['1 and 2';'1 and 4'])
        end
        cnt=cnt+1;
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([units(units_cnt).name(1:end-27), '    Crosscorrelations'],'Interpreter','None','FontSize',15,'FontWeight','bold')
    
    saveas(gcf,['F:\',date,'\processing\crosscorr\Corr_',units(units_cnt).name(1:end-27),'.bmp'])
  
end

% convolved firing rate
length(units)
for sig=1:3:15
    
    for units_cnt=1:10
        
        load([unitspath,units(units_cnt).name])
        
        
        cnt=1;
        close (figure(1))
        figure(1)
        set(gcf,'Position', [1 31 1680 946])
        for i=1:4:32
            flips=spike_info.flip_times{i}(1:1800:end,1); % flips in ms, rounded
            trial1=i;
            trial2=i+1;
            tmp=spike_info.spike_times{trial1};
            spikes1=0;
            for j=2:2:length(flips)
                spikes1=spikes1+(convolved(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999,sig,30000))/4;
            end
            spikes1=spikes1(121:end-120);
            tmp=spike_info.spike_times{trial2};
            spikes2=0;
            for j=2:2:length(flips)
                spikes2=spikes2+(convolved(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999,sig,30000))/4;
            end
            spikes2=spikes2(121:end-120);
            [c_ww,lags] = xcorr(spikes1,spikes2,500,'coef');
            subplot(2,4,cnt)
            plot(lags,c_ww,'linewidth',2)
            hold on
            
            
            
            
            trial1=i;
            trial2=i+3;
            tmp=spike_info.spike_times{trial1};
            spikes1=0;
            for j=2:2:length(flips)
                spikes1=spikes1+(convolved(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999,sig,30000))/4;
            end
            spikes1=spikes1(121:end-120);
            tmp=spike_info.spike_times{trial2};
            spikes2=0;
            for j=2:2:length(flips)
                spikes2=spikes2+(convolved(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999,sig,30000))/4;
            end
            spikes2=spikes2(121:end-120);
            
            [c_ww,lags] = xcorr(spikes1,spikes2,500,'coef');
            subplot(2,4,cnt)
            plot(lags,c_ww,'r','linewidth',2)
            
            
            
            
            trial1=i+2;
            trial2=i+3;
            tmp=spike_info.spike_times{trial1};
            spikes1=0;
            for j=2:2:length(flips)
                spikes1=spikes1+(convolved(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999,sig,30000))/4;
            end
            spikes1=spikes1(121:end-120);
            tmp=spike_info.spike_times{trial2};
            spikes2=0;
            for j=2:2:length(flips)
                spikes2=spikes2+(convolved(tmp(tmp<flips(j)&tmp>flips(j-1)+1000)-flips(j-1)-999,sig,30000))/4;
            end
            spikes2=spikes2(121:end-120);
            
            [c_ww,lags] = xcorr(spikes1,spikes2,500,'coef');
            subplot(2,4,cnt)
            plot(lags,c_ww,'color',[0 0.7 0],'linewidth',2)
            
            
            
            
            title(['ND', int2str(9-cnt)])
            if cnt==1
                legend(['1 and 2';'1 and 4';'3 and 4'])
            end
            line([0,0],get(gca,'YLim'),'color','k','linewidth',2)
            cnt=cnt+1;
        end
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title([units(units_cnt).name(1:end-27), '    Crosscorrelations of convolved spike rate (SIG=',int2str(sig),') averaged within each trial (only high contrast)'],'Interpreter','None','FontSize',15,'FontWeight','bold')
        if ~exist(['F:\',date,'\processing\crosscorr_sig',int2str(sig)])
            mkdir(['F:\',date,'\processing\crosscorr_sig',int2str(sig)]);
        end
        saveas(gcf,['F:\',date,'\processing\crosscorr_sig',int2str(sig),'\Corr_',units(units_cnt).name(1:end-27),'.bmp'])
        
    end
end

%raster
date='20120329';
unitspath=['F:\',date,'\processing\easy_formatted_units\'];
units=dir([unitspath,'*FFFlicker.mat']);
for units_cnt=1:length(units)
    
    load([unitspath,units(units_cnt).name])
    close (figure(1))
    figure(1)
    a=subplot(1,1,1);
    set(gcf,'Position', [1 31 1680 946])
    cnt=0;tmp1=[];
    for i=5:32
        tmp=spike_info.spike_times{i};
        flips=spike_info.flip_times{i}([3600 5400]);
        tmp=tmp(tmp<flips(2)&tmp>flips(1))-flips(1)-5000;
        tmp(tmp<0)=[];
        tmp(tmp>20000)=[];
        length(tmp)
        tmp1=[tmp1 tmp+cnt*20000];
        cnt=cnt+1;
    end
    rasterplot(tmp1,28,20000,a)
    for i=1:6
        line([0,30000],[5.75+(i-1)*6,5.75+(i-1)*6],'color','m')
    end
    set(gca,'YTick',3:6:40,'YTickLabel',{'ND7','ND6','ND5','ND4','ND3','ND2','ND1'})
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([units(units_cnt).name(1:end-27), '    2nd repetition of high contrast, 5s-25s'],'Interpreter','None','FontSize',15,'FontWeight','bold')
    
    saveas(gcf,['F:\',date,'\processing\rasters\Raster_',units(units_cnt).name(1:end-27),'.bmp'])
  
end



%cross binarized
date='20120329';
unitspath=['F:\',date,'\processing\easy_formatted_units\'];
units=dir([unitspath,'*FFFlicker.mat']);
bin=20;
for units_cnt=1:length(units)
    
    load([unitspath,units(units_cnt).name])
    spikes=zeros(32,floor(31000/bin));
    for i=1:32
        flips=spike_info.flip_times{i}(1:1800:end,1); % flips in ms, rounded
        tmp=spike_info.spike_times{i};
        spikes1=zeros(1,31000);
        for j=2:2:length(flips)
            spikes1(tmp(tmp<flips(j)&tmp>flips(j-1))-flips(j-1))=spikes1(tmp(tmp<flips(j)&tmp>flips(j-1))-flips(j-1))+1;
        end
        spikes1=reshape(spikes1(1:floor(31000/bin)*bin),bin,31000/bin);
        spikes(i,:)=sum(spikes1(:,1:floor(31000/bin)));        
    end
    a=sum(spikes);    
    spikes(:,a==0)=[];
    
    close (figure(1))
    figure(1)
    set(gcf,'Position', [1 31 1680 946])
    cnt=1;
    for i=1:4:32
        subplot(2,4,(i+3)/4)
        hold on
        [c_ww,lags] = xcorr(spikes(i+2,:),spikes(i+3,:),50,'coef');
        plot(lags,c_ww,'color',[0 0.7 0],'linewidth',2)
        [c_ww,lags] = xcorr(spikes(i,:),spikes(i+2,:),50,'coef');
        plot(lags,c_ww,'linewidth',2)
        [c_ww,lags] = xcorr(spikes(i,:),spikes(i+3,:),50,'coef');
        plot(lags,c_ww,'r','linewidth',2)        
        title(['ND', int2str(9-cnt)])
        if cnt==1
            legend(['1 and 2';'1 and 4';'3 and 4'])
        end
        line([0,0],get(gca,'YLim'),'color','k')
        cnt=cnt+1;
    end
    
%     HeatMap(corr(spikes'))
        
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([units(units_cnt).name(1:end-27), '    Crosscorrelation, high contrast mean, bin ',int2str(bin),'ms'],'Interpreter','None','FontSize',15,'FontWeight','bold')
    
    saveas(gcf,['F:\',date,'\processing\crosscorr_bin20\Crosscorr_',units(units_cnt).name(1:end-27),'.bmp'])
  
end




%cross binarized positions of maxima
date='20120329';
unitspath=['F:\',date,'\processing\easy_formatted_units\'];
units=dir([unitspath,'*FFFlicker.mat']);
bin=20;
for units_cnt=1:length(units)
    units_cnt
    load([unitspath,units(units_cnt).name])
    spikes=zeros(32,floor(31000/bin));
    for i=1:32
        flips=spike_info.flip_times{i}(1:1800:end,1); % flips in ms, rounded
        tmp=spike_info.spike_times{i};
        spikes1=zeros(1,31000);
        for j=2:2:length(flips)
            spikes1(tmp(tmp<flips(j)&tmp>flips(j-1))-flips(j-1))=spikes1(tmp(tmp<flips(j)&tmp>flips(j-1))-flips(j-1))+1;
        end
        spikes1=reshape(spikes1(1:floor(31000/bin)*bin),bin,31000/bin);
        spikes(i,:)=sum(spikes1(:,1:floor(31000/bin)));        
    end
    a=sum(spikes);    
    spikes(:,a==0)=[];
    
    cnt=1;
    for i=1:4:32
        [c_ww,lags] = xcorr(spikes(i+2,:),spikes(i+3,:),50,'coef');
        [~,b]=max(c_ww);
        max_pos(units_cnt,cnt,1)=lags(b);
        [c_ww,lags] = xcorr(spikes(i,:),spikes(i+2,:),50,'coef');
        [~,b]=max(c_ww);
        max_pos(units_cnt,cnt,2)=lags(b);
        [c_ww,lags] = xcorr(spikes(i,:),spikes(i+3,:),50,'coef');
        [~,b]=max(c_ww);
        max_pos(units_cnt,cnt,3)=lags(b);
        cnt=cnt+1;
    end

end


close (figure(1))
figure(1)
set(gcf,'Position', [1 31 1680 946])
for i=1:8
    subplot(2,4,i)
    hold on
    a=reshape(max_pos(:,i,:),52,3);
    plot(a(:,1),'b*')
    plot(60,mean(a(abs(a(:,1))<5,1)),'b.','Markersize',10)
    plot(a(:,2),'ro')
    plot(65,mean(a(abs(a(:,2))<5,2)),'r.','Markersize',10)
    plot(a(:,3),'gx')
    plot(70,mean(a(abs(a(:,3))<5,3)),'g.','Markersize',10)
    axis([0 75 -5 5])
    line([0,75],[0,0],'color','k')
end

