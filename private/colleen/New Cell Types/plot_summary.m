%% ------- INPUTS -------
date='2008-12-12-1';
concatname='data006-nwpca/data006';
staruns='data006'; % long WN runs for STA low to high NDF

cell_specification = [{'OFF parasol'}, {'OFF large 3'}];

%% ------ END OF INPUTS --------
% file path to save pictures
filepath=['/Users/colleen/Desktop/Summary/',date,'/',concatname,'/'];
if ~exist(filepath,'dir')
    mkdir(filepath);
end

% load STA runs
starun = cell(1,size(staruns,1));
for i=1:size(staruns,1)
    fullname=[staruns(i,:)];
    datarun = load_data(fullfile(server_path(),date,concatname,fullname));
    datarun = load_params(datarun,'verbose',1);
    datarun = load_neurons(datarun);
    datarun = load_sta(datarun);
    datarun = set_polarities(datarun);
    starun = datarun;
end



fig=figure('PaperPosition',[0 0 10 9],'PaperSize',[10 9]);
set(fig,'color','white','position',[82 242 1785 856]);
% 0.05 0.55, 0.6,0.4
% make plots
for i=1:size(cell_specification,2)
    [cell_numbers] = get_cell_indices(starun, cell_specification{i});
    
    % Vision ID of the cell
    visionID = cell_specification{i};
%     ind = cell_numbers(i);
    
    %
    
    % plot stuff
    
    
    % plot mosaic
    h=subplot('position',[0.05+ 0.3*(i-1), 0.75, 0.20,0.20]);
    %     rasterplot(nmrasters,(num_repeats-1) * size(nmruns,1),movie_length,h)
    plot_axes = plot_rf_summaries(datarun, cell_specification{i})
    axis equal
    if i == 1
        title([cell_specification{1},' Mosaic'])
    else
        title([cell_specification{2},' Mosaic'])
    end
    
    hold on
    h=subplot('position',[0.05+0.3*(i-1), 0.45 0.20,0.20]);
    sta_tc=zeros(size(starun.vision.timecourses(2).r,1),size(starun,1));
    sta_tc_bins=zeros(size(starun.vision.timecourses(2).r,1),size(starun,1));
    for k = 1:size(cell_numbers,2)
        if ~isempty(starun.vision.timecourses(cell_numbers(k)).g)
            sta_tc(:,k*3-2)=starun.vision.timecourses(cell_numbers(k)).r;
            sta_tc(:,k*3-1)=starun.vision.timecourses(cell_numbers(k)).g;
            sta_tc(:,k*3)=starun.vision.timecourses(cell_numbers(k)).b;
                    
            sta_tc(:,(k*3-2):k*3) =sta_tc(:,(k*3-2):k*3)./norm(sta_tc(:,(k*3-2):k*3));
            
            tc_bin=starun.stimulus.refresh_period;
            
            
        end
    end
    sta_tc_bins(:,1)=-tc_bin*(size(starun.vision.timecourses(2).r,1)-3):tc_bin:tc_bin*2;
    red = sta_tc(:,1:3:end);
    green = sta_tc(:,2:3:end);
    blue = sta_tc(:,3:3:end);
    plot(h,sta_tc_bins,red,'Color', [1 0 0 0.15]);
    hold on
    plot(h,sta_tc_bins,blue,  'Color', [0 0 1 0.15]);
    plot(h,sta_tc_bins,green, 'Color', [0 1 0 0.15]);
    plot(h, sta_tc_bins, zeros(size(sta_tc_bins)), 'k')
    axis tight
    title('Timecourses')
end

sta_coords = [0.05 0.27, 0.12,0.12;...
    0.05 0.14, 0.12,0.12;...
    0.05 0.01, 0.12,0.12;...
    0.2 0.27, 0.12,0.12;...
    0.2 0.14, 0.12,0.12;...
    0.2 0.01, 0.12,0.12;...
    0.35 0.27, 0.12,0.12;...
    0.35 0.14, 0.12,0.12;...
    0.35 0.01, 0.12,0.12;...
    0.50 0.27, 0.12,0.12;...
    0.50 0.14, 0.12,0.12;...
    0.50 0.01, 0.12,0.12;...
    0.65 0.27, 0.12,0.12;...
    0.65 0.14, 0.12,0.12;...
    0.65 0.01, 0.12,0.12;...
    0.80 0.27, 0.12,0.12;...
    0.80 0.14, 0.12,0.12;...
    0.80 0.01, 0.12,0.12];

if size(cell_numbers,2) > 18
sta_coords = [0.05 0.27, 0.10,0.10;...
    0.05 0.14, 0.10,0.10;...
    0.05 0.01, 0.10,0.10;...
    0.16 0.27, 0.10,0.10;...
    0.16 0.14, 0.10,0.10;...
    0.16 0.01, 0.10,0.10;...
    0.27 0.27, 0.10,0.10;...
    0.27 0.14, 0.10,0.10;...
    0.27 0.01, 0.10,0.10;...
    0.38 0.27, 0.10,0.10;...
    0.38 0.14, 0.10,0.10;...
    0.38 0.01, 0.10,0.10;...
    0.49 0.27, 0.10,0.10;...
    0.49 0.14, 0.10,0.10;...
    0.49 0.01, 0.10,0.10;...
    0.60 0.27, 0.10,0.10;...
    0.60 0.14, 0.10,0.10;...
    0.60 0.01, 0.10,0.10;...
    0.71 0.27, 0.10,0.10;...
    0.71 0.14, 0.10,0.10;...
    0.71 0.01, 0.10,0.10;...
    0.82 0.27, 0.10,0.10;...
    0.82 0.14, 0.10,0.10;...
    0.82 0.01, 0.10,0.10;...
    0.82 0.66, 0.10,0.10;...
    0.82 0.53, 0.10,0.10;...
    0.82 0.40, 0.10,0.10;...
    0.71 0.40, 0.10,0.10;...
    0.71 0.53, 0.10,0.10;...
    0.71 0.66, 0.10,0.10];
end

% plot fine STA
for k=1:size(cell_numbers,2)
    h=subplot('position',sta_coords(k,:));
    sta=squeeze(datarun.stas.stas{cell_numbers(k)});
    colormap gray
    if size(sta, 3) ~= 3
        [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1));
        
        imagesc(sta(:,:,start_index))
    else
        [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));
        
        sta = norm_image(sta);
        image(sta(:,:,:,start_index))
    end
    
    set(h,'xtick',0,'ytick',0)
%      tmp = size(sta);
%         axis([0.5 tmp(1)+0.5 0.5 tmp(2)+0.5])
end





% save figure
print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,'/'],[cell_specification{2}]));
%     close(fig)



