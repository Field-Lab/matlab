
date='2015-03-09-2';
concatname='d05-27-norefit';
finalname='-from-d05-d27';
runs=['data006'; 'data009'; 'data012'; 'data016'; 'data020'; 'data024'];
staruns=['data005'; 'data008'; 'data014'; 'data018'; 'data022'; 'data026'];
% NDF5
% data005 WN
% data006 NSEM
% data007 WN repeats
% 
% NDF4 
% data008 WN
% data009 NSEM
% data010 WN repeats
% 
% NDF3
% data014 WN
% data012 NSEM
% data013 WN repeats
% 
% NDF2
% data018 WN
% data016 NSEM
% data017 WN repeats
% 
% NDF1
% data022 WN
% data020 NSEM
% data021 WN repeats
% 
% NDF0
% data026 WN
% data024 NSEM
% data025 WN repeats

filepath=['/Users/alexth/Desktop/Light_adaptation/movie_rasters_fr_comb/',date,'/',concatname,'/'];
if ~exist(filepath,'dir')
    mkdir(filepath);
end

% NSEM runs
runpaths=[];
for i=1:size(runs,1)
    fullname=[runs(i,:),finalname];
    runpaths=[runpaths; fullfile(server_path(),date,concatname,fullname,fullname)];
end
  
datarun = load_data(runpaths(1,:));
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

datarun2 = load_data(runpaths(2,:));
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_neurons(datarun2);

datarun3 = load_data(runpaths(3,:));
datarun3 = load_params(datarun3,'verbose',1);
datarun3 = load_neurons(datarun3);

datarun4 = load_data(runpaths(4,:));
datarun4 = load_params(datarun4,'verbose',1);
datarun4 = load_neurons(datarun4);

datarun5 = load_data(runpaths(5,:));
datarun5 = load_params(datarun5,'verbose',1);
datarun5 = load_neurons(datarun5);

datarun6 = load_data(runpaths(6,:));
datarun6 = load_params(datarun6,'verbose',1);
datarun6 = load_neurons(datarun6);

%STA runs
runpaths=[];
for i=1:size(staruns,1)
    fullname=[staruns(i,:),finalname];
    runpaths=[runpaths; fullfile(server_path(),date,concatname,fullname,fullname)];
end

starun = load_data(runpaths(1,:));
starun = load_params(starun,'verbose',1);
starun = load_sta(starun);
starun = set_polarities(starun);

starun2 = load_data(runpaths(2,:));
starun2 = load_params(starun2,'verbose',1);
starun2 = load_sta(starun2);
starun2 = set_polarities(starun2);

starun3 = load_data(runpaths(3,:));
starun3 = load_params(starun3,'verbose',1);
starun3 = load_sta(starun3);
starun3 = set_polarities(starun3);

starun4 = load_data(runpaths(4,:));
starun4 = load_params(starun4,'verbose',1);
starun4 = load_sta(starun4);
starun4 = set_polarities(starun4);

starun5 = load_data(runpaths(5,:));
starun5 = load_params(starun5,'verbose',1);
starun5 = load_sta(starun5);
starun5 = set_polarities(starun5);

starun6 = load_data(runpaths(6,:));
starun6 = load_params(starun6,'verbose',1);
starun6 = load_sta(starun6);
starun6 = set_polarities(starun6);


ncells=length(datarun.cell_ids);

trigs=datarun.triggers;
figure
plot(diff(trigs))
myTrigs=find(diff(trigs)>0.9&diff(trigs)<2);
myTrigs=myTrigs(myTrigs<740);
myTrigs=[0; myTrigs];



for i=1:ncells
    
    for j=1:length(starun.cell_types)
        if ~isempty(find(datarun.cell_ids(i)==starun.cell_types{j}.cell_ids, 1))
            folder=starun.cell_types{j}.name;
            if ~exist([filepath,folder],'dir')
                mkdir([filepath,folder]);
            end
            break
        end
    end
    
    spikes{1}=datarun.spikes{i};
    spikes{2}=datarun2.spikes{i};
    spikes{3}=datarun3.spikes{i};
    spikes{4}=datarun4.spikes{i};
    spikes{5}=datarun5.spikes{i};
    spikes{6}=datarun6.spikes{i};
    
    splitSpikes=cell(19,6);
    
    for j=1:19
        for k=1:6
            tmp=spikes{k}(spikes{k}>trigs(myTrigs(j)+1) & spikes{k}<trigs(myTrigs(j+1)))...
                - trigs(myTrigs(j)+1);
            splitSpikes{j,k}=tmp*1000;
        end
    end
    
    t=[];
    cnt=[0 20 40 60 80 100];
    fr=zeros(6,31000);
    splitRasters=[];
    for j=1:19
        for k=1:6
            splitRasters=[splitRasters splitSpikes{j,k}'+31000*cnt(k)];             
            tmp=convolved(splitSpikes{j,k}',40,31000);
            fr(k,:)=fr(k,:)+tmp(((size(tmp,2)-31000)/2+1):end-((size(tmp,2)-31000)/2))/20;        
        end
        cnt=cnt+1;
    end
    splitRasters=sort(splitRasters);

    fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[82 242 1785 856]);
    
    
    h=subplot('position',[0.05 0.55, 0.6,0.4]);
    rasterplot(splitRasters,cnt(6),31000,h)
    for k=1:5
        line([0,31000], [28.25,28.25]+30*(k-1),'color','r','linewidth',0.7)
    end
    for k=1:29
        line([0,0]+k*1000,[0,178],'color','b','linewidth',0.3)
    end
    set(gca,'ytick',14:30:178,'yticklabel',{'NDF5.0','NDF4.0','NDF3.0','NDF2.0','NDF1.0','NDF0.0'})
    set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'})
    axis([0 30000 0 Inf])
    ylabel('light level')
    xlabel('NSEM,s')
    title(['2015-03-09-2, cell ',int2str(datarun.cell_ids(i))])
    
    h=subplot('position',[0.05 0.05, 0.6,0.4]);
    hold on
    plot(fr','linewidth',0.7)
    axis([0 30000 0 Inf])
    set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'})
    legend('5.0','4.0','3.0','2.0','1.0','0.0');
    
    
    h=subplot('position',[0.67 0.73, 0.12,0.23]);
    sta=squeeze(starun.stas.stas{i});
    colormap gray
    imagesc(sta(:,:,26))
    plot_rf_summaries(starun, starun.cell_ids(i), 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    title('NDF 5.0, 26')
    
    h(2)=subplot('position',[0.67 0.48, 0.12,0.23]);
    sta=squeeze(starun2.stas.stas{i});
    colormap gray
    imagesc(sta(:,:,26))
    hold on
    plot_rf_summaries(starun2, starun2.cell_ids(i), 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    title('NDF 4.0, 26')
    
    h(3)=subplot('position',[0.67 0.23, 0.12,0.23]);
    sta=squeeze(starun3.stas.stas{i});
    colormap gray
    imagesc(sta(:,:,26))
    plot_rf_summaries(starun3, starun3.cell_ids(i), 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    title('NDF 3.0, 26')
    
    h(4)=subplot('position',[0.84 0.73, 0.12,0.23]);    
    sta=squeeze(starun4.stas.stas{i});
    colormap gray
    imagesc(sta(:,:,26)) 
    plot_rf_summaries(starun4, starun4.cell_ids(i), 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    title('NDF 2.0, 26')
    
    h(5)=subplot('position',[0.84 0.48, 0.12,0.23]);
    sta=squeeze(starun5.stas.stas{i});
    colormap gray
    imagesc(sta(:,:,28))
    plot_rf_summaries(starun5, starun5.cell_ids(i), 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    title('NDF 1.0, 28')
    
    h(6)=subplot('position',[0.84 0.23, 0.12,0.23]);
    sta=squeeze(starun6.stas.stas{i});
    colormap gray
    imagesc(sta(:,:,27)) 
    plot_rf_summaries(starun6, starun6.cell_ids(i), 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    title('NDF 0.0, 27')
    
    for j=1:6
        set(h(j),'xtick',0,'ytick',0)
%         axis([0.5 16.5 0.5 16.5])
    end

    
    h=subplot('position',[0.67 0.05, 0.32,0.17]);
    sta_tc=zeros(30,6);
    sta_tc_bins=zeros(30,6);
    if ~isempty(starun.vision.timecourses(i).g)
        sta_tc(:,1)=starun.vision.timecourses(i).g;
        tc_bin=starun.stimulus.refresh_period;
        sta_tc_bins(:,1)=-tc_bin*27:tc_bin:tc_bin*2;
    end
    if ~isempty(starun2.vision.timecourses(i).g)
        sta_tc(:,2)=starun2.vision.timecourses(i).g;
        tc_bin=starun2.stimulus.refresh_period;
        sta_tc_bins(:,2)=-tc_bin*27:tc_bin:tc_bin*2;
    end
    if ~isempty(starun3.vision.timecourses(i).g)
        sta_tc(:,3)=starun3.vision.timecourses(i).g;
        tc_bin=starun3.stimulus.refresh_period;
        sta_tc_bins(:,3)=-tc_bin*27:tc_bin:tc_bin*2;
    end
    if ~isempty(starun4.vision.timecourses(i).g)
        sta_tc(:,4)=starun4.vision.timecourses(i).g;
        tc_bin=starun4.stimulus.refresh_period;
        sta_tc_bins(:,4)=-tc_bin*27:tc_bin:tc_bin*2;
    end
    if ~isempty(starun5.vision.timecourses(i).g)
        sta_tc(:,5)=starun5.vision.timecourses(i).g;
        tc_bin=starun5.stimulus.refresh_period;
        sta_tc_bins(:,5)=-tc_bin*27:tc_bin:tc_bin*2;
    end
    if ~isempty(starun6.vision.timecourses(i).g)
        sta_tc(:,6)=starun6.vision.timecourses(i).g;
        tc_bin=starun6.stimulus.refresh_period;
        sta_tc_bins(:,6)=-tc_bin*27:tc_bin:tc_bin*2;
    end

    plot(sta_tc_bins,sta_tc, 'linewidth',2);
    legend('5.0','4.0','3.0','2.0','1.0','0.0','location','northwest')
    hold on
    line([min(sta_tc_bins(:)) max(sta_tc_bins(:))],[0,0],'color','k')
    axis tight
   
        
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(datarun.cell_ids(i))]));    
    close(fig)

end

