
date='2015-02-24-5';
% concatname='d02_07_13_16-norefit';
% runs=['data002'; 'data007'; 'data013'; 'data016']; 
concatname='d01_02_03_06_07_08_12_13_14_15_16_17-norefit';
finalname='-from-data001_data002_data003_data006_data007_data008_data012_data013_data014_data015_data016_data017';
runs=['data002'; 'data007'; 'data013'; 'data016'];
staruns=['data001'; 'data006'; 'data012'; 'data015'];

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
    
    splitSpikes=cell(19,4);
    
    for j=1:19
        for k=1:4
            tmp=spikes{k}(spikes{k}>trigs(myTrigs(j)+1) & spikes{k}<trigs(myTrigs(j+1)))...
                - trigs(myTrigs(j)+1);
            splitSpikes{j,k}=tmp*1000;
        end
    end
    
    t=[];
    cnt=[0 20 40 60];
    fr=zeros(4,31000);
    splitRasters=[];
    for j=1:19
        for k=1:4
            splitRasters=[splitRasters splitSpikes{j,k}'+31000*cnt(k)];             
            tmp=convolved(splitSpikes{j,k}',40,31000);
            fr(k,:)=fr(k,:)+tmp(((size(tmp,2)-31000)/2+1):end-((size(tmp,2)-31000)/2))/20;        
        end
        cnt=cnt+1;
    end
    splitRasters=sort(splitRasters);

    figure
    set(gcf,'position',[82 242 1785 856])
    h=subplot('position',[0.05 0.55, 0.6,0.4]);
    rasterplot(splitRasters,cnt(4),31000,h)
    for k=1:3
        line([0,31000], [28.25,28.25]+31*(k-1),'color','r','linewidth',1.2)
    end
    for k=1:29
        line([0,0]+k*1000,[0,117],'color','b','linewidth',0.8)
    end
    set(gca,'ytick',14:28:98,'yticklabel',{'NDF3.0','NDF2.0'})
    set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'})
    axis([0 30000 0 Inf])
    ylabel('light level')
    xlabel('NSEM,s')
    title(['2015-02-24-5, cell ',int2str(datarun.cell_ids(i))])
    
    h=subplot('position',[0.05 0.05, 0.6,0.4]);
    hold on
    plot(fr','linewidth',1.2)
    axis([0 30000 0 Inf])
    set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'})
    legend('NDF 3.0', 'NDF 2.0','NDF 3.0', 'NDF 2.0');
    
    
    h=subplot('position',[0.67 0.63, 0.15,0.32]);
    sta=squeeze(starun.stas.stas{i});
    colormap gray
    imagesc(sta(:,:,28))
    plot_rf_summaries(starun, starun.cell_ids(i), 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    title('NDF 3.0 first')
    
    h(2)=subplot('position',[0.67 0.26, 0.15,0.32]);
    sta=squeeze(starun2.stas.stas{i});
    colormap gray
    imagesc(sta(:,:,28))
    hold on
    plot_rf_summaries(starun2, starun2.cell_ids(i), 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    title('NDF 2.0 first')
    
    h(3)=subplot('position',[0.84 0.63, 0.15,0.32]);
    sta=squeeze(starun3.stas.stas{i});
    colormap gray
    imagesc(sta(:,:,28))
    plot_rf_summaries(starun3, starun3.cell_ids(i), 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    title('NDF 3.0 second')
    
    h(4)=subplot('position',[0.84 0.26, 0.15,0.32]);    
    sta=squeeze(starun4.stas.stas{i});
    colormap gray
    imagesc(sta(:,:,28)) 
    plot_rf_summaries(starun4, starun4.cell_ids(i), 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    title('NDF 2.0 second')
    
    for j=1:4
        set(h(j),'xtick',0,'ytick',0)
        axis([0.5 16.5 0.5 16.5])
    end
    
    h=subplot('position',[0.67 0.05, 0.32,0.17]);
    sta_tc=zeros(30,4);
    if ~isempty(starun.vision.timecourses(i).g)
        sta_tc(:,1)=starun.vision.timecourses(i).g;
    end
    if ~isempty(starun2.vision.timecourses(i).g)
        sta_tc(:,2)=starun2.vision.timecourses(i).g;
    end
    if ~isempty(starun3.vision.timecourses(i).g)
        sta_tc(:,3)=starun3.vision.timecourses(i).g;
    end
    if ~isempty(starun4.vision.timecourses(i).g)
        sta_tc(:,4)=starun4.vision.timecourses(i).g;
    end
    plot(sta_tc);
    legend('3.0 first','2.0 first','3.0 second','2.0 second','location','northwest')
    hold on
    line([0 30],[0,0],'color','k')
    tmp=starun.stimulus.refresh_period;
    ticks=sort(30:-200/tmp:0);
    set(h,'xtick',ticks,'xticklabel',{'-1800','-1600','-1400','-1200','-1000',...
        '-800','-600','-400','-200','0'})
    
    saveas(gcf,[filepath,folder,'/cell_',int2str(datarun.cell_ids(i)),'.jpg'])
    close(gcf)
end

