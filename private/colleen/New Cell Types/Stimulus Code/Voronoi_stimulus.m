%% --------------- Voronoi_stimulus (focal white noise) -------------------
%% Function: Generate a mask that can be copied to the stimulus machines to generate white noise just over a particular region (or particular regions) to target individual cells

%% How to use: 
%     1) In vision, select one or more cells that you want white noise to be place over (the rest will be gray). It is best if they aren't touching
%     2) Enter the file name of the white noise classification run the RF will be based on
%     3) Check the file path to save the mask to (will change when using computers in the lab)
%     4) screen size should stay at 640x480
%     5) Enter the stixel size, which has to match the white noise run as well as the stixel size of the Voroni stimulus

%% Potential problems:
%     1) No significant stixels are found
%         -mess with significant stixel code parameters
%     2) I think the stimulus has to be the same size stixels as the classification run the mapping comes from
    

%% Inputs
% cells : cell IDs from Vision
% file_name : start with piece number and continue until the .sta etc files (eg '2006-06-06-2/data003/data003')
% mdf_file : white noise xml for this classification run
% file_path : where you want the mask saved
% screen_size_y : height of mask (keep at 480)
% screen_size_x : width of mask (keep at 640)
% stixel_size : match the stixel size of the classification run

%% Results
% Graph of Vision STA for the selected cells and their time courses
% Graph of what the mask will look like. Colored pixels will have white
% noise
% Mask in a text file saved to the indicated location

%% Author 
% Colleen Rhoades (rhoades@stanford.edu)
% April 7, 2015


clear
%% ------------------------------ INPUTS -----------------------------------
cells = {5329};
file_name = '2015-04-09-5/data001/data001';
mdf_file='/Volumes/Analysis-1/stimuli/white-noise-xml/BW-20-8-0.48-22222-16x16.xml';
file_path = ['/Users/vision/matlab/private/colleen/New Cell Types/Stimulus Code/2014-04-09-5/data001_', num2str(cells{1}), '.txt'];
screen_size_y =320; % vertical size
screen_size_x = 320; % hortizontal size
stixel_size = 20;
%% ------------------------------- Load Data ------------------------------------------

datarun.names.rrs_params_path=['/Volumes/Acquisition/Analysis/', file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Acquisition/Analysis/', file_name, '.sta'];
opt=struct('verbose',1,'load_params',1,'load_neurons',0,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',false);
opt.load_sta_params.frames = 1:30;% if this line is missing, will error; have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

%% ------------------------------- Plot Vision STA -----------------------------
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
    %% ------------------------- Find RF location -------------------------
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
%% --------------------------- Compute size of white noise patch (matches for all cells) ---------------------
% Chose largest width and height
% Add to the right edge for cells missing odd number of stixels
diff_x = overall_size(:,2)-overall_size(:,1)+1;
diff_y = overall_size(:,4)-overall_size(:,3)+1;
loc_x = max(diff_x);
loc_y = max(diff_y);
for k = 1:size(overall_size,1)
    current = overall_size(k,:);
    mean_x = mean([current(1), current(2)]);
    mean_y = mean([current(3), current(4)]);
    x_stix(k,:) = round(mean_x-(loc_x-1)/2:mean_x+(loc_x-1)/2);
    y_stix(k,:) = round(mean_y-(loc_y-1)/2:mean_y+(loc_y-1)/2);
end

%% ---------------------------- Get all important stixels  --------------------------
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

%% --------------------------------- Check if cells are on the screen ------------------------
if max(large_cell_stixels_final(:,1)*stixel_size) > screen_size_x
    disp('error: x dimension of chosen stixels doesn''t fit on screen')
    return
elseif max(large_cell_stixels_final(:,2)*stixel_size) > screen_size_y
    disp('error: y dimension of chosen stixels doesn''t fit on screen')
    return
end


%% -------------------------------- Populate the mask --------------------------------------
scale_factor_x = screen_size_x/ stixel_size;
scale_factor_y = screen_size_y/stixel_size;
subdivide = 2;
count = 1;
for i = 1:size(large_cell_stixels_final,1)
    
    subdivide_factor = mod(i,4);
    pix_y =     stixel_size*large_cell_stixels_final(i, 1) - stixel_size/subdivide + 1: stixel_size*large_cell_stixels_final(i, 1);
    pix_x =    stixel_size*large_cell_stixels_final(i, 2) - stixel_size/subdivide + 1: stixel_size*large_cell_stixels_final(i, 2);
    
    myMap(pix_x,pix_y) = count;
    count = count+1;
    pix_y =     -stixel_size/subdivide+ stixel_size*large_cell_stixels_final(i, 1) - stixel_size/subdivide + 1: stixel_size*large_cell_stixels_final(i, 1)- stixel_size/subdivide;
    pix_x =    -stixel_size/subdivide+ stixel_size*large_cell_stixels_final(i, 2) - stixel_size/subdivide + 1: stixel_size*large_cell_stixels_final(i, 2) - stixel_size/subdivide;
    
    myMap(pix_x,pix_y) = count;
    count = count+1;
    
    pix_y =     stixel_size*large_cell_stixels_final(i, 1) - stixel_size/subdivide + 1: stixel_size*large_cell_stixels_final(i, 1);
    pix_x =    -stixel_size/subdivide+ stixel_size*large_cell_stixels_final(i, 2) - stixel_size/subdivide + 1: stixel_size*large_cell_stixels_final(i, 2) - stixel_size/subdivide;
    myMap(pix_x,pix_y) = count;
    count = count+1;
    
    pix_y =     -stixel_size/subdivide+ stixel_size*large_cell_stixels_final(i, 1) - stixel_size/subdivide + 1: stixel_size*large_cell_stixels_final(i, 1)- stixel_size/subdivide;
    pix_x =    stixel_size*large_cell_stixels_final(i, 2) - stixel_size/subdivide + 1: stixel_size*large_cell_stixels_final(i, 2);
    myMap(pix_x,pix_y) = count;
    count = count+1;
    
end



%% ------------------------------- Write the mask to a file ---------------------------------
dlmwrite(file_path, myMap, 'delimiter', '\t', 'newline', 'pc'); % if this errors which save path
savedMap = dlmread(file_path);

%% ------------------------------- Display mask ----------------------------------------------
figure
imagesc(savedMap)
title('Stixels to be modulated')
axis equal