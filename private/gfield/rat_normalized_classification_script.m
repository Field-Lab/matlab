%% This script is for loading in rat data from matlab and playing around
% with classifying the RFC types

%% load datarun
datarun = load_data('/snle/lab/Experiments/Array/Analysis/2011-09-14-1/data000/data000');
datarun = load_params(datarun);
datarun = load_sta(datarun);

% compute the time courses and RFs for every cell
datarun = get_sta_summaries(datarun, 'all');

% identify a cell type of interest
cell_type = 'OFF type 1';

% get the cell indices for this cell type
cell_indices = get_cell_indices(datarun, cell_type);

% get how many cells are in this list
num_cells = length(cell_indices);

% initialize matrix of time course
num_frames = length(datarun.stas.time_courses{1});
time_courses = zeros(num_cells, num_frames);

% loop through this list of cells and put time courses in time_courses
for rgc = 1:num_cells
    
    time_courses(rgc,:) = datarun.stas.time_courses{cell_indices(rgc)};
    
end


% calculate and plot the mean time course for the cells in this list
mean_time_course = mean(time_courses);

figure(1)
plot(mean_time_course)


%% Now we want to normalize this time course to peak at 1

monitor_refresh = 120; % in units of Hz
time_bin_size = datarun.stimulus.interval * (1/monitor_refresh) * 1000; % units of ms
total_tc_length = -1*(num_frames -1)*time_bin_size;

time_points = total_tc_length:time_bin_size:0


% make a nice fig showing the time course with labeled and meaningful axes
figure(1)
plot(time_points, mean_time_course)
xlabel('ms')
ylabel('relative contrast')
title('mean time course')

% get the peak frame
[~, peak_frame_index] = max(abs(mean_time_course));

normalized_time_points = time_points ./ abs(time_points(peak_frame_index));

% make a fig showing this normalized
figure(2)
plot(normalized_time_points, mean_time_course)
xlabel('ms')
ylabel('relative contrast')
title('mean time course')
