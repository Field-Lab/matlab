%% Identify data of interest
% set main path to data
main_path = '/Volumes/lab/Experiments/Array/Analysis/'; % location in directory to get data
master_path_tail = '2011-06-30-0/data003/data003';
slave_path_tail = '2011-06-30-0/data004/data004';

slave_path_tail = '2011-06-30-0/data005/data005';

%% Load Data

% define cell types of interest
cell_types = {1,2,3,4};

% data, piece, datarun of MASTER white noise run
master_full_path = [main_path, master_path_tail];
% load data from MASTER
opt=struct('verbose',1,'load_params',1,'load_neurons',1, 'load_ei',1);
datarun{1} = load_data(master_full_path, opt);
datarun{1} = load_sta(datarun{1}, 'load_sta', 'all', 'save_sta', false);
datarun{1} = get_sta_summaries(datarun{1}, cell_types, 'keep_stas', false, 'keep_rfs', false);

% path to slave white noise repeat data
slave_full_path = [main_path, slave_path_tail];
datarun{2} = load_data(slave_full_path, opt);
% map cells
%datarun=map_cell_types(datarun); % for classical mapping
[cell_list_map, failed_cells] = map_ei(datarun{1}, datarun{2}, 'master_cell_type', cell_types,...
                                'slave_cell_type', 'all', 'corr_threshold', 0.95);

%NOTE: cell_list_map lists every cell_id in slave that corresponds to the index in Master.
% for example, if the first cell_id in MASTER corresponds to cell_ID 31 in
% slave, then the first listed number in cell_list_map will be 31.
                            
%% make raster plot
% get relevant trigger set 

master_ind = get_cell_indices(datarun{1}, datarun{1}.cell_types{3}.cell_ids(2));
slave_ind = get_cell_indices(datarun{2}, cell_list_map{master_ind});


temp_cell_type = {1,2,3,4};
clear slave_cell_ids temp_indices master_cell_ids
temp_indices = get_cell_indices(datarun{1}, temp_cell_type);
mapped_cell_indices = [];
for rgc = 1:length(temp_indices);
    if isempty(cell_list_map{temp_indices(rgc)})
        slave_cell_ids(rgc) = 0;
    else
        slave_cell_ids(rgc) = cell_list_map{temp_indices(rgc)};
        mapped_cell_indices = [mapped_cell_indices, rgc];
    end
end
slave_cell_ids = slave_cell_ids(mapped_cell_indices);
master_cell_ids = datarun{1}.cell_ids(temp_indices(mapped_cell_indices));
datarun{2} = compute_rasters(datarun{2}, slave_cell_ids,...
            'triggers_per_trial', 6, 'master_cell_ids', master_cell_ids,...
            'print_summaries', true, 'save_path', '/Analysis/raster-files/all-cells/',...
            'bin_size', 0.01);
        

%% Load Cone files and compute generator for WN repeats        
        
cone_file_path = '/Volumes/lab/Experiments/Array/Shared/one/2011-06-30-0_rf-3-bayes-msf_200.00-BW-2-4/';
datarun{1} = import_single_cone_data(datarun{1},cone_file_path);
cd(cone_file_path)
load Wc

% put white noise movie into datarun{1} because it has "right" triggers
wn_repeat_movie_path = '/Volumes/lab/acquisition/movie-xml/BW-2-4-0.48-11111-300x300-60.35.xml';
datarun{1} = load_java_movie(datarun{1},wn_repeat_movie_path);

% now get just the frames that were shown on each trial of the white noise
% repeats
stimulus_interval = 4; %frames
stimulus_duration = 594; % frames... just under 10 seconds worth of frames
white_noise_frames = floor(stimulus_duration/stimulus_interval);
cone_inputs = zeros(white_noise_frames,size(Wc,2));
field_width = 300;
field_height = 300;

% project each frame through the cone receptive fields (in Wc)
fprintf('Computing cone input in frames %d to %d... \n',0,white_noise_frames)
for ss = 1:white_noise_frames

    % get new frame
    STAFrame = datarun{1}.stimulus.java_movie.getFrame(ss-1);
    new_frame = permute(reshape(STAFrame.getBuffer,3,field_width,field_height),[3 2 1]) - .5;

    new_frame = reshape(new_frame,[],1);

    % convert to cone space
    cone_inputs(ss,:) = full(Wc'*double(new_frame));

end

datarun{1}.mapping.map_to_data004 = cell_list_map;
datarun{2}.cones.cone_inputs = cone_inputs;


%% remap to a new datarun

mapped_datarun.spikes = datarun{2}.spikes;

        








