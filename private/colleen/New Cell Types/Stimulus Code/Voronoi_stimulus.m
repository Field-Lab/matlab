
clear

% Code to make stimulus to just have white noise over the large cells and a
% gray background everywhere else


%% ------------------------------ INPUTS -----------------------------------
cells = 2689;
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
opt.load_sta_params.frames = 1:30% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);
[cell_numbers, cell_type, cell_type_number] = get_cell_indices(datarun, cells);
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
num_stixels = length(num_x)*length(num_y);
large_cell_stixels = zeros(num_stixels,2);
n = length(num_y);
for i =1 :length(num_x)
    large_cell_stixels(i*n-n+1:n*i, 1) = repmat(num_x(i), n,1);
end
large_cell_stixels(:,2) = repmat(num_y',length(num_x),1);
large_cell_stixels = fliplr(large_cell_stixels);
% large_cell_stixels =[y,x]; % horizontal over from top left then vertical down from top corner

% large_cell_stixels = [10,5; 10,10; 2,12]; % horizontal over from top left then vertical down from top corner
if max(large_cell_stixels(:,1)*stixel_size) > screen_size_x
    disp('error: x dimension of chosen stixels doesn''t fit on screen')
    return
elseif max(large_cell_stixels(:,2)*stixel_size) > screen_size_y
    disp('error: y dimension of chosen stixels doesn''t fit on screen')
    return
end

    

myMap = zeros(screen_size_y, screen_size_x); % pixesl on the screen
scale_factor_x = screen_size_x/ stixel_size;
scale_factor_y = screen_size_y/stixel_size;

for i = 1:size(large_cell_stixels,1)
    pix_y = (stixel_size*large_cell_stixels(i, 1)-(stixel_size-1)):(stixel_size)*large_cell_stixels(i,1);
    pix_x = (stixel_size*large_cell_stixels(i,2)-(stixel_size-1)):(stixel_size)*large_cell_stixels(i,2);
    myMap(pix_x,pix_y) = i;
end



dlmwrite(file_path, myMap, 'delimiter', '\t', 'newline', 'pc');
savedMap = dlmread(file_path);
figure
imagesc(savedMap)
title('Stixels to be modulated')
axis equal