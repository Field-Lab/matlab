% Plot large cells across NDFs
clear all
close all 
clc

date='2008-04-08-2';
concatname='d00_03_05_07';
staruns={'data000'; 'data003';'data005'; 'data007'}; % long WN runs for STA
final_name = '-from-data000_data003_data005_data007';
% wnruns=['data006'; 'data019';]; % WN repeats
cell_specification = [652 947 1277 1775 1896 2423 5944 6067 7490];
% cell_specification = {'ON large 1'};

type = {'ON large 2'};%cell_specification;
NDFs_wn = {'4.3', '4', '0.6', '0.6'};
titles = {'NDF 4.3', 'NDF 4', 'NDF 0.6', 'NDF 0.6'};

filepath=['/Users/colleen/Desktop/Light_adaptation/NDFs',date,'/',concatname,'/'];
if ~exist(filepath,'dir')
    mkdir(filepath);
end

% load STA runs
starun = cell(1,size(staruns,1));
for i=1:size(staruns,1)
    fullname=cell2mat([staruns(i,:), final_name]);
    datarun = load_data(fullfile(server_path(),date,concatname,fullname,fullname));
    datarun = load_params(datarun,'verbose',1);
    datarun = load_sta(datarun);
    datarun = set_polarities(datarun);
    starun{i} = datarun;
end

% 
% sta_coords = [0.67 0.73, 0.12,0.23;...
%     0.67 0.48, 0.12,0.23;...
%     0.67 0.23, 0.12,0.23;...
%     0.84 0.73, 0.12,0.23;...
%     0.84 0.48, 0.12,0.23;...
%     0.84 0.23, 0.12,0.23];

[cell_numbers] = get_cell_indices(starun{1}, cell_specification);
vision_ids = datarun.cell_ids(cell_numbers);

    
    
 % plot stuff
    fig=figure('PaperPosition',[0 0 12*length(cell_numbers)/10 8],'PaperSize',[12*length(cell_numbers)/10 8]);
    set(fig,'color','white','position',[82 242 1785*length(cell_numbers)/10 856]);
    
 


for i=1:length(cell_numbers)
   
    % Vision ID of the cell
    visionID = vision_ids(i);
    ind = cell_numbers(i);
    % identify cell type and create folder

    if i == 1
        if iscellstr(cell_specification)
            folder = cell_specification{1};
        else
            [folder, ~] = find_cell_type(starun{1}, cell_specification(i));
        end
        if ~exist([filepath,folder],'dir')
            mkdir([filepath,folder]);
        end
    end
    
    
    
   
    % plot fine STA
    for k=1:size(staruns,1)
        hold on
%         h=subplot('position',sta_coords(k,:));
        h = subplot(size(NDFs_wn,2), length(vision_ids), sub2ind([length(vision_ids), size(NDFs_wn,2)],i,k));
        sta=squeeze(starun{k}.stas.stas{ind});
        colormap gray
        if size(sta, 3) ~= 3
            [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1));
            
            imagesc(sta(:,:,start_index))
        else
            [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));
            
            sta = norm_image(sta);
            image(sta(:,:,:,start_index))
        end
        
        plot_rf_summaries(starun{k}, visionID, 'clear', false,  'plot_fits', true, 'fit_color', 'r')
%         title('NDF 5.0, 26')
        set(h,'xtick',0,'ytick',0)
        tmp = size(sta);
        axis([0.5 tmp(2)+0.5 0.5 tmp(1)+0.5])
        if mod(k,size(NDFs_wn,2)) == 1
%             xlabel(['Cell ', num2str(vision_ids(i))])
            title({['Cell ', num2str(vision_ids(i))];titles{k}}, 'fontsize', 10)

        else
                    title(titles{k}, 'fontsize', 10)

        end
    end
   
    
%     % plot time course from vision
%     h=subplot('position',[0.7 0.03, 0.29,0.15]);
%     sta_tc=zeros(30,5);
%     sta_tc_bins=zeros(30,5);
%     for k = 1:size(staruns,1)
%         if ~isempty(starun{k}.vision.timecourses(ind).g)
%             sta_tc(:,k)=starun{k}.vision.timecourses(ind).g;
%             tc_bin=starun{k}.stimulus.refresh_period;
%             sta_tc_bins(:,k)=-tc_bin*27:tc_bin:tc_bin*2;
%         end
%     end
%     plot(sta_tc_bins,sta_tc, 'linewidth',2);
%     legend(NDFs_wn,'location','northwest')
%     hold on
%     line([min(sta_tc_bins(:)) max(sta_tc_bins(:))],[0,0],'color','k')
%         axis tight
%     current_xlim = get(gca, 'xlim');
% 
%     set(gca, 'xlim', [-1000, current_xlim(2)])
%     title('Coarse Timecourses')
    
    % save figure
%     print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(visionID)]));

    
end
%         close(fig)
suptitle({date; type{1}})
hgexport(gcf, sprintf('%s%s%s.pdf',[filepath,folder,'/'],[type{1}]))

% export_fig(sprintf('%s%s%s.pdf',[filepath,folder,'/'],[type{1}]), '-pdf')
