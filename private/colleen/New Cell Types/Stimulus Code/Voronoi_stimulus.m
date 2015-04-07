
clear

% Code to make stimulus to just have white noise over the large cells and a
% gray background everywhere else


%% ------------------------------ INPUTS -----------------------------------
cells = {2689, 4067};
file_name = '2006-06-06-2/data012/data012';
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111.xml';
file_path = '/Users/colleen/matlab/private/colleen/New Cell Types/Stimulus Code/test.m';
screen_size_y =320; % vertical size
screen_size_x = 640; % hortizontal size
stixel_size = 10;
%%


% cell_specification = [502,860,1024,1130,2076,2361,2618,2705,3022,3172,3213,3559,4022,4071,4238,4774,4852,5496,6518,6533,6860,7279,7671];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];

datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];
opt=struct('verbose',1,'load_params',1,'load_neurons',0,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',false);
% opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames = 1:30;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);
myMap = zeros(screen_size_y, screen_size_x); % pixesl on the screen

for cell = 1:size(cells,2)
    [cell_numbers, cell_type, cell_type_number] = get_cell_indices(datarun, cells{cell});
    sta = datarun.stas.stas{cell_numbers};
    sig_stixels = significant_stixels(sta);
    % sig_stixels = significant_stixels(sta, 'select', 'thresh', 'thresh', 4)
    
    % find peak frame
    time_course = time_course_from_sta(sta, sig_stixels);
    
    %% find the peak and trough of TC
    
    % collapse TC to one color if RGB
    if size(time_course, 2) > 1;
        summed_tc = sum(time_course,2); % collapse RGB tc
    else
        summed_tc = time_course; % don't collapse a BW tc
    end
    
    % matlab indexes the time course such that the highest index is the frame
    % closest to the spike time. The fit function assumes the opposite.
    % To compensate the min(max) frame is subtracted from the total frame
    % number
    frame_number = size(sta,4);
    [max_tc_val, max_frame] = max(summed_tc);
    [min_tc_val, min_frame] = min(summed_tc);
    max_frame = frame_number - max_frame;
    min_frame = frame_number - min_frame;
    
    % Sort parameters appropriately for ON and OFF cells
    if min_frame < max_frame % true for off cells
        peak_frame = min_frame;
    else                     % true for on cells
        peak_frame = max_frame;
    end
    figure;
    if size(time_course, 2) > 1;
        imagesc(squeeze(sta(:,:,2,frame_number - peak_frame)))
    else
        imagesc(squeeze(sta(:,:,1,frame_number - peak_frame)))
    end
    
    title('Peak STA from Vision (Green)')
    axis equal
    
    %%
    
    [x,y] = find(full(sig_stixels == 1));
    min_x = min(x);
    min_y = min(y);
    max_x = max(x);
    max_y = max(y);
    num_x = min_x:max_x;
    num_y = min_y:max_y;
    overall_size(cell, :) = [min_x, max_x, min_y, max_y];
    num_stixels(cell,:) = length(num_x)*length(num_y);
    
end
diff_x = overall_size(:,2)-overall_size(:,1)+1;
diff_y = overall_size(:,4)-overall_size(:,3)+1;
loc_x = max(diff_x);
loc_y = max(diff_y);
for k = 1:size(overall_size,1)
    current = overall_size(k,:)
    mean_x = mean([current(1), current(2)])
    mean_y = mean([current(3), current(4)])
    x_stix(k,:) = round(mean_x-(loc_x-1)/2:mean_x+(loc_x-1)/2);
    y_stix(k,:) = round(mean_y-(loc_y-1)/2:mean_y+(loc_y-1)/2);




end
        large_cell_stixels = zeros(length(x_stix(1,:)).*length(y_stix(1,:)),2);

for k = 1:size(overall_size,1)
    x_stix2 = x_stix(k,:);
    y_stix2 = y_stix(k,:);
    n = length(y_stix2);
    for i =1 :length(x_stix2)
        large_cell_stixels(i*n-n+1:n*i, 1) = repmat(x_stix2(i), n,1);
    end
    large_cell_stixels(1*[1:loc_x*loc_y],2) = repmat(y_stix2',length(x_stix2),1);
    large_cell_stixels_final(k*[loc_x*loc_y]- [loc_x*loc_y] + 1:k*[loc_x*loc_y], :) = large_cell_stixels(1*[1:loc_x*loc_y], :);
end

    large_cell_stixels_final = fliplr(large_cell_stixels_final);

    % large_cell_stixels =[y,x]; % horizontal over from top left then vertical down from top corner
    
    % large_cell_stixels = [10,5; 10,10; 2,12]; % horizontal over from top left then vertical down from top corner
    if max(large_cell_stixels_final(:,1)*stixel_size) > screen_size_x
        disp('error: x dimension of chosen stixels doesn''t fit on screen')
        return
    elseif max(large_cell_stixels_final(:,2)*stixel_size) > screen_size_y
        disp('error: y dimension of chosen stixels doesn''t fit on screen')
        return
    end
    
    
    
    scale_factor_x = screen_size_x/ stixel_size;
    scale_factor_y = screen_size_y/stixel_size;
    for i = 1:size(large_cell_stixels_final,1)
       
        pix_y = (stixel_size*large_cell_stixels_final(i, 1)-(stixel_size-1)):(stixel_size)*large_cell_stixels_final(i,1);
        pix_x = (stixel_size*large_cell_stixels_final(i,2)-(stixel_size-1)):(stixel_size)*large_cell_stixels_final(i,2);
        myMap(pix_x,pix_y) = i;
    end
    



dlmwrite(file_path, myMap, 'delimiter', '\t', 'newline', 'pc');
savedMap = dlmread(file_path);
figure
imagesc(savedMap)
title('Stixels to be modulated')
axis equal