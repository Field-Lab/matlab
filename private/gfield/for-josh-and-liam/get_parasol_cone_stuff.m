%% LOAD DATA
datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/Chichilnisky-lab/2008-08-27-5/data003/data003/data003');
                fruit_name = 'plantain/';
                num_frames = 144103;
                num_cones = 1870;
                num_rgcs = 376;
                tc_length = 7;
                data_path = '/Volumes/lab/Experiments/1cone/plantain/';
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/plantain/'];     

%%                
% load movie_xml_path & other info
datarun = load_index(datarun, 'index_path', '/Volumes/lab/Experiments/Array/Analysis/Chichilnisky-lab/2008-08-27-5/Index');

% load spikes times, trigger times, params, stas and cone info
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, datarun.names.nickname);
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

datarun = get_sta_summaries(datarun, {1,2}, 'keep_stas', false);

%%

% load cone weights matrix
load([data_path, '/Wc.mat'])
  
% --- read cone generator signals from disk ---
disp('Reading file cone_input.bin...')
file_name = [data_path, 'cone_input.bin'];
fid = fopen(file_name, 'r');
cone_inputs = fread(fid, [num_frames, num_cones], 'float32', 0,'b');
fclose(fid);

% --- read spike_rate file from disk ---
disp('Reading file spike_count.bin...')
file_name = [data_path, 'spike_counts.bin'];
fid = fopen(file_name, 'r');
spike_rate = fread(fid, [num_frames, num_rgcs], 'float32', 0,'b');
spike_rate = spike_rate';
fclose(fid);

% --- read the time courses for the RGCs ---
disp('Reading time course file...')
file_name = [data_path, 'rgc_tcs.bin'];
time_course_fid = fopen(file_name, 'r');
rgc_tcs = fread(time_course_fid, [tc_length, num_rgcs],'float32', 0, 'b');
rgc_tcs = rgc_tcs';
fclose(time_course_fid);

% -- load RGC cell IDs ---
disp('reading cell IDs...')
file_name = [data_path,'rgc_IDs.txt'];
rgc_ids = dlmread(file_name, '\t');

% --- load cone IDs ---
disp('reading cone IDs...')
file_name = [data_path, 'cone_IDs.txt'];
cone_ids = dlmread(file_name, '\t');

% go through each cell and make spike_time cell from spike_rate matrix
for cc = 1:length(rgc_ids)

    % translate to spike times (with duplicates for multiple spikes per time bin)
    spike_times_ = [];
    for nn = 1:max(spike_rate(cc,:))
        spike_times_ = [spike_times_, find( spike_rate(cc,:) > (nn-1) )];
    end
    
    % put into storage variables
    spike_times{cc} = sort(spike_times_);
    
end


datarun.cones.Wc = Wc;
datarun.cones.cone_inputs = cone_inputs;
datarun.cones.spike_rate = spike_rate;
datarun.cones.spike_times = spike_times;
datarun.cones.rgc_tcs = rgc_tcs;
datarun.cones.rgc_ids = rgc_ids;
datarun.cones.cone_ids = cone_ids;


