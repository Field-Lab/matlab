
date='2015-03-09-2';
concatname='d05-27-norefit';
finalname='-from-d05-d27';
nmruns=['data006'; 'data009'; 'data012'; 'data016'; 'data020'; 'data024']; % NSEM runs
staruns=['data005'; 'data008'; 'data014'; 'data018'; 'data022'; 'data026']; % long WN runs for STA
wnruns=['data007'; 'data010'; 'data013'; 'data017'; 'data021'; 'data025']; % WN repeats

% file path to save pictures
filepath=['/Users/alexth/Desktop/Light_adaptation/NSEM_WN/',date,'/',concatname,'/'];
if ~exist(filepath,'dir')
    mkdir(filepath);
end

% load NSEM runs
nmrun = cell(1,6);
for i=1:size(nmruns,1)
    fullname=[nmruns(i,:),finalname];
    datarun = load_data(fullfile(server_path(),date,concatname,fullname,fullname));
    datarun = load_params(datarun,'verbose',1);
    datarun = load_neurons(datarun);
    nmrun{i} = datarun;
end

% load STA runs
starun = cell(1,6);
for i=1:size(staruns,1)
    fullname=[staruns(i,:),finalname];
    datarun = load_data(fullfile(server_path(),date,concatname,fullname,fullname));
    datarun = load_params(datarun,'verbose',1);
    datarun = load_sta(datarun);
    datarun = set_polarities(datarun);
    starun{i} = datarun;
end

% load WN repeats runs
wnrun = cell(1,6);
for i=1:size(wnruns,1)
    fullname=[wnruns(i,:),finalname];
    datarun = load_data(fullfile(server_path(),date,concatname,fullname,fullname));
    datarun = load_params(datarun,'verbose',1);
    datarun = load_neurons(datarun);
    wnrun{i} = datarun;
end

clear datarun fullname finalname

% triggers NSEM - most likely the same, but just in case
nm_trigs = zeros(20,6);
for i=1:6
    nm_trigs(:,i)=[0; find(diff(nmrun{i}.triggers)>0.9)];
end

% triggers NSEM - most likely the same, but just in case
wn_trigs = zeros(20,6);
for i=1:6
    wn_trigs(:,i)=[0; find(diff(wnrun{i}.triggers)>0.84)];
end

% STA subplots coordinates and relevant frames
% instead frame, one could make a convolution of the STA with the time course
sta_coords = [0.67 0.73, 0.12,0.23;...
    0.67 0.48, 0.12,0.23;...
    0.67 0.23, 0.12,0.23;...
    0.84 0.73, 0.12,0.23;...
    0.84 0.48, 0.12,0.23;...
    0.84 0.23, 0.12,0.23];
sta_frames = [26 26 26 26 28 27];


% make plots
for i=1:length(nmrun{1}.cell_ids)
    
    % Vision ID of the cell
    visionID = starun{1}.cell_ids(i);
    
    % identify cell type and create folder  
    [folder, ~] = find_cell_type(starun{1}, visionID);
    if ~exist([filepath,folder],'dir')
        mkdir([filepath,folder]);
    end
    
    
    % NSEM processing    
    nmrasters = [];
    cnt = 0; % 1st trial number for each NDF
    movie_length = 30000; % in ms
    for k=1:6
        spikes = nmrun{k}.spikes{i};
        trigs = nmrun{k}.triggers;
        beg_points = nm_trigs(:,k);
        for j=1:19
            tmp=spikes(spikes>trigs(beg_points(j)+1) & spikes<trigs(beg_points(j+1)))...
                - trigs(beg_points(j)+1);
            nmrasters=[nmrasters tmp'*1000 + movie_length*cnt];
            cnt = cnt+1;
        end
    end
    
    % WN repeats processing
    wnrasters = [];
    cnt = 0; % 1st trial number for each NDF
    wn_length = 30000; % in ms
    for k=1:6
        spikes = wnrun{k}.spikes{i};
        trigs = wnrun{k}.triggers;
        beg_points = wn_trigs(:,k);
        for j=1:19
            tmp=spikes(spikes>trigs(beg_points(j)+1) & spikes<trigs(beg_points(j+1)))...
                - trigs(beg_points(j)+1);
            wnrasters=[wnrasters tmp'*1000 + wn_length*cnt];
            cnt = cnt+1;
        end
    end
    
    
    % plot stuff
    fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[82 242 1785 856]);    
    
    % plot NSEM
    h=subplot('position',[0.05 0.55, 0.6,0.4]);
    rasterplot(nmrasters,19 * 6,movie_length,h)
    for k=1:5
        line([0,30000], [19,19]*k*1.5,'color','r','linewidth',0.7)
    end
    for k=1:29
        line([0,0]+k*1000,[0,19*6*1.5],'color','b','linewidth',0.3)
    end
    set(gca,'ytick',19*0.75:19*1.5:19*1.5*6,'yticklabel',{'5','4','3','2','1','0'})
    set(gca,'xtick',5000:5000:25000,'xticklabel',{'5','10','15','20','25'})
    axis([0 30000 0 Inf])
    ylabel('light level, NDF')
    xlabel('NSEM,s')
    title(['2015-03-09-2, cell ',int2str(visionID)])
    
    
    % plot WN repeats
    h=subplot('position',[0.05 0.09, 0.6,0.4]);
    rasterplot(wnrasters,19 * 6,wn_length,h)
    for k=1:5
        line([0,30000], [19,19]*k*1.5,'color','r','linewidth',0.7)
    end
    set(gca,'ytick',19*0.75:19*1.5:19*1.5*6,'yticklabel',{'5','4','3','2','1','0'})
    set(gca,'xtick',5000:5000:25000,'xticklabel',{'5','10','15','20','25'})
    axis([0 30000 0 Inf])
    ylabel('light level, NDF')
    xlabel('WN repeats,s')
    
    % plot STA
    for k=1:6
        h=subplot('position',sta_coords(k,:));
        sta=squeeze(starun{k}.stas.stas{i});
        colormap gray
        imagesc(sta(:,:,sta_frames(k)))
        plot_rf_summaries(starun{k}, visionID, 'clear', false,  'plot_fits', true, 'fit_color', 'r')
        title('NDF 5.0, 26')
        set(h,'xtick',0,'ytick',0)
        tmp = size(sta);
        axis([0.5 tmp(1)+0.5 0.5 tmp(2)+0.5])
        title(['NDF',int2str(6-k),', frame ',int2str(sta_frames(k))])
    end
    
  
    % plot time course from vision      
    h=subplot('position',[0.67 0.05, 0.32,0.17]);
    sta_tc=zeros(30,6);
    sta_tc_bins=zeros(30,6);    
    for k = 1:6
        if ~isempty(starun{k}.vision.timecourses(i).g)
            sta_tc(:,k)=starun{k}.vision.timecourses(i).g;
            tc_bin=starun{k}.stimulus.refresh_period;
            sta_tc_bins(:,k)=-tc_bin*27:tc_bin:tc_bin*2;
        end
    end
    plot(sta_tc_bins,sta_tc, 'linewidth',2);
    legend('5','4','3','2','1','0','location','northwest')
    hold on
    line([min(sta_tc_bins(:)) max(sta_tc_bins(:))],[0,0],'color','k')
    axis tight   
    
    % save figure
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(visionID)]));    
    close(fig)

end

