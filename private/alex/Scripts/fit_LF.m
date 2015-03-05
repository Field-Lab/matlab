




% fitting
% parameters: a1 - amplitude, b1 - delay (position of the peak), c1 - width
clear
date='20120920';
filter_length=500;
filterspath=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/FFFlicker_LF/'];
filters=dir([filterspath,'*filter.mat']);

load([filterspath,filters(1).name])
trials=size(HighFilters,2);

on=0;off=0;
intervals=100;

lower_border=-10000;
upper_border=10000;
spikeLowerBorder=50;
peaks_inds=zeros(length(filters),trials);
a1=zeros(length(filters),trials);
b1=zeros(length(filters),trials);
c1=zeros(length(filters),trials);
on_off=zeros(length(filters),1)-10;
names=cell(1,length(filters));



[ampl_max, ind_max]=max(HighFilters(1:200,:));
[ampl_min, ind_min]=min(HighFilters(1:200,:));




for filters_cnt=1:length(filters)
    filters_cnt
    load([filterspath,filters(filters_cnt).name])
    names{filters_cnt}=filters(filters_cnt).name;
    [~, ind_max]=max(HighFilters(1:200,:));
    [~, ind_min]=min(HighFilters(1:200,:));
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

save(['F:\',date,'\processing\linear_filters_fit_v1\gaussian_parameters.mat'],'a1','b1','c1','names','peaks_inds')

% plotting fit
load(['F:\',date,'\processing\linear_filters_fit_v1\gaussian_parameters.mat'])
for filters_cnt=1:length(filters)
    if abs(on_off(filters_cnt))==1
        load([filterspath,filters(filters_cnt).name])
        close (figure(1))
        figure(1)
        set(gcf,'Position', [1 31 1680 946])
        for trial=1:100
            subplot(10,10,trial)
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

figure
for i=1:12:12*8
    plot(i:i+11,b1(3,i:i+11),'*')
    hold on
    line([i+11.5,i+11.5],[0, max(b1(3,:))])
end
figure
for i=1:12:12*8
    plot(i:i+11,a1(3,i:i+11),'*')
    hold on
    line([i+11.5,i+11.5],[min(a1(3,:)),0])
end

