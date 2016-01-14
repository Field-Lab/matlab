% Plot large cells across NDFs
clear all
close all
clc

date='2015-09-23-7';
concatname='d19-39';
staruns={'data019'; 'data029';'data030'; 'data031';'data032'; 'data035'; 'data037'; 'data038'; 'data020'; 'data023';'data025';'data026';'data028';'data039'}; % long WN runs for STA
final_name = '-from-data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039';
% wnruns=['data006'; 'data019';]; % WN repeats
cell_specification = [530 1103 2378 2644 4043 5299 6788 918 3931 4132 5075 5506 6830 7203 532 1006 1010 1597 1880 2421 3398 4658 5342 5493 7087 7145 81 1175 3394 4657 4850 4864 5090 5361 6323 7131];
% cell_specification = {'ON large 1'};

type = {'OFF large 1'};%cell_specification;
NDFs_wn = {'0', '0', '0', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2'};
titles = {'NDF 0 8-1', 'NDF 0 8-8', 'NDF 0 4-2', 'NDF 2 4-4', 'NDF 2 4-12', 'NDF 2 4-4','NDF 2 4-4','NDF 2 4-4','NDF 2 4-4','NDF 2 4-4','NDF 2 4-4','NDF 2 4-4','NDF 2 4-4','NDF 2 4-4'};

filepath=['/Users/colleen/Desktop/Light_adaptation/NDFs',date,'/',concatname,'/'];
% if ~exist(filepath,'dir')
%     mkdir(filepath);
% end

% load STA runs
starun = cell(1,size(staruns,1));
for i=1:size(staruns,1)
    fullname=cell2mat([staruns(i,:), final_name]);
    datarun = load_data(fullfile(server_path(),date,concatname,fullname,fullname), 'load_neurons', 1);
    datarun = load_params(datarun,'verbose',1);
    datarun = load_sta(datarun);
    datarun = load_spikes(datarun);

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
        spike_count(k,i) = length(starun{k}.spikes(ind))
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
