%% ------- INPUTS -------
date='2005-04-26-0';
concatname='data009/';
staruns='data009'; % long WN runs for STA low to high NDF
padding= 10;
cell_specification = [{'OFF parasol' }, {'OFF parasol'}];
num_frames = 20;
vision =0;
%% ------ END OF INPUTS --------
% file path to save pictures
filepath=['/Users/colleen/Desktop/GrantJuly2015/',date,'/',concatname,'/'];
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
% make plots
for i=1:size(cell_specification,2)
    [cell_numbers] = get_cell_indices(starun, cell_specification{i});
    
    % Vision ID of the cell
    visionID = cell_specification{i};
    %     ind = cell_numbers(i);
    
    %
    
    % plot stuff
    
    
    % plot mosaic
    if vision == 1
        h=subplot('position',[0.05+ 0.3*(i-1), 0.75, 0.20,0.20]);
        plot_axes = plot_rf_summaries(datarun, cell_specification{i});
    end
    
    
    
    if vision == 1
    sta_tc=zeros(60,size(starun,1));
    else
            sta_tc=zeros(num_frames,size(starun,1));

    end
    
    sta_tc_bins=zeros(num_frames,size(starun,1));
    for k = 1:size(cell_numbers,2)
        
        %% from fitting code
        if vision~=1
            fit_params = datarun2.matlab.sta_fits{cell_numbers(k)}.fit_params;
            fixed_params = datarun2.matlab.sta_fits{cell_numbers(k)}.fixed_params;
            fit_indices = datarun2.matlab.sta_fits{cell_numbers(k)}.fit_indices;
            fixed_indices = datarun2.matlab.sta_fits{cell_numbers(k)}.fixed_indices;
            temp_stix = datarun2.matlab.sta_fits{cell_numbers(k)}.sig_stixels;
            
            all_params(fit_indices) = fit_params;
            all_params(fixed_indices) = fixed_params;
            sta_fit = sta_fit_function(all_params);
            
            %
            h=subplot('position',[0.05+ 0.3*(i-1), 0.75, 0.20,0.20]);
            hold on
            h(k) = plot_spatial_sd(all_params);
            hold off
            axis equal
        end
        
        if i == 1
            title([cell_specification{1},' Mosaic'])
        else
            title([cell_specification{2},' Mosaic'])
        end
        if ~isempty(starun.vision.timecourses(cell_numbers(k)).g)
            if vision == 1
                sta_tc(:,k*3-2)=starun.vision.timecourses(cell_numbers(k)).r;
                sta_tc(:,k*3-1)=starun.vision.timecourses(cell_numbers(k)).g;
                sta_tc(:,k*3)=starun.vision.timecourses(cell_numbers(k)).b;
                
                sta_tc(:,(k*3-2):k*3) =sta_tc(:,(k*3-2):k*3)./norm(sta_tc(:,(k*3-2):k*3));
            end
            
            tc_bin=starun.stimulus.refresh_period;
            
            
            if vision~=1
                fit_tc_one = time_course_from_sta(sta_fit, temp_stix);
                sta_tc(:,k*3-2) = fit_tc_one(length(fit_tc_one)-num_frames+1:length(fit_tc_one),1);
                sta_tc(:,k*3-1) = fit_tc_one(length(fit_tc_one)-num_frames+1:length(fit_tc_one),2);
                sta_tc(:,k*3) = fit_tc_one(length(fit_tc_one)-num_frames+1:length(fit_tc_one),3);
                %             norm_factor = max(abs(reshape(fit_tc, 1, [])));
                %             fit_tc = fit_tc ./ norm_factor;
                sta_tc(:,(k*3-2):k*3) =sta_tc(:,(k*3-2):k*3)./norm(sta_tc(:,(k*3-2):k*3));
            end
            
            
        end
    end
    sta_tc = sta_tc(size(sta_tc,1)-num_frames+1:size(sta_tc,1),:);
    
    sta_tc_bins(:,1)=-tc_bin*(num_frames -2):tc_bin:tc_bin;
    red = sta_tc(:,1:3:end);
    green = sta_tc(:,2:3:end);
    blue = sta_tc(:,3:3:end);
    hold on
    h=subplot('position',[0.05+0.3*(i-1), 0.45 0.20,0.20]);
    plot(h,sta_tc_bins,red,'Color', [1 0 0 0.15]);
    hold on
    plot(h,sta_tc_bins,blue,  'Color', [0 0 1 0.15]);
    plot(h,sta_tc_bins,green, 'Color', [0 1 0 0.15]);
    plot(h, sta_tc_bins, zeros(size(sta_tc_bins)), 'k')
    axis tight
    xlabel('Time to spike (ms)')
    ylabel('STA (arb. units)')
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
    axis equal
    if vision~=1
        bounds = autozoom_to_fit(datarun2, cell_numbers(k), padding, [1,1], 1, vision);
    else
        bounds = autozoom_to_fit(datarun, cell_numbers(k), padding, [1,1], 1, vision);
    end
    
    sta=squeeze(datarun.stas.stas{cell_numbers(k)});
    colormap gray
    if size(sta, 3) ~= 3
        [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1));
        
        image(sta(:,:,start_index))
        
    else
        [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));
        
        sta = norm_image(sta);
        image(sta(:,:,:,start_index))
    end
    
    if vision ~= 1
        the_fit.mean  = [datarun2.matlab.sta_fits{cell_numbers(k)}.center_point_x  datarun2.matlab.sta_fits{cell_numbers(k)}.center_point_y];
        the_fit.sd  = [datarun2.matlab.sta_fits{cell_numbers(k)}.center_sd_x  datarun2.matlab.sta_fits{cell_numbers(k)}.center_sd_y];
        the_fit.angle = datarun2.matlab.sta_fits{cell_numbers(k)}.center_rotation_angle;
    else
        the_fit = datarun.stas.fits{cell_numbers(k)};
    end
    
    % get points of an ellipse with these parameters
    [X,Y] = drawEllipse([the_fit.mean the_fit.sd the_fit.angle]);
    hold on
    plot(X,Y,'Color','k', 'LineWidth',1)
    % %
    
    
    axis([bounds])
    axis square
    axis off
    pix_field_width = datarun.stimulus.field_width*datarun.stimulus.stixel_width;
    scale_bar_width = 15; %pixels
    stix_scale_bar = (scale_bar_width/pix_field_width)*datarun.stimulus.field_width;
    plot([bounds(1)+1; bounds(1)+stix_scale_bar+1], [bounds(4)-1; bounds(4)-1], '-k', 'LineWidth', 2)
    
    %     axis([[10  34]    -60  140])
    
    
    %     set(h,'xtick',0,'ytick',0)
    %      tmp = size(sta);
    %         axis([0.5 tmp(1)+0.5 0.5 tmp(2)+0.5])
end





% save figure
print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,'/'],[cell_specification{2}]));
%     close(fig)



