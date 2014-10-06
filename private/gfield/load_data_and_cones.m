function [datarun, cone_info] = load_data_and_cones(data_name, varargin)

p = inputParser;
p.addRequired('data_name', @ischar);
p.addParamValue('sta_summaries', false)

p.parse('data_name', varargin{:});



%% LOAD DATA
switch data_name

    case 'apple';  datarun = load_data('2010-03-05-2', 'rf-13-apple-gf');
                fruit_name = 'apple/';
                num_frames = 76284;
                num_cones = 2287;
                num_rgcs = 295;
                tc_length = 6;
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/apple/'];
                
    case 'peach';  datarun = load_data('2008-08-27-0','rf-1-peach');
                fruit_name = 'peach/';
                num_frames = 48015;
                num_cones = 2107;
                num_rgcs = 328;
                tc_length = 7;
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/peach/'];                

    case 'plantain';  datarun = load_data('2008-08-27-5', 'rf-3-plantain');
                fruit_name = 'plantain/';
                num_frames = 144103;
                num_cones = 1870;
                num_rgcs = 376;
                tc_length = 7;
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/plantain/'];     
    case 'blueberry';  datarun = load_data('2008-08-26-2','rf-1-blueberry');
                fruit_name = 'blueberry/';
                data_path = [single_cone_path,fruit_name];
                num_frames = 63046;
                num_cones = 1382;
                num_rgcs = 226;
                tc_length = 6;
                midget_window_size = 22;
                parasol_window_size = 50;
                save_path = ['~/Desktop/blueberry/'];
    case 'apricot';  datarun = load_data('2009-04-13-5', 'rf-5-apricot');
                fruit_name = 'apricot/';
                data_path = [single_cone_path,fruit_name];
                num_frames = 108078;
                num_cones = 6197;
                num_rgcs = 837;
                tc_length = 7;
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/apricot/'];
    case 'kiwi';  datarun = load_data('2008-05-13-3','rf-6-kiwi');
                fruit_name = 'kiwi/';
                data_path = [single_cone_path,fruit_name];
                num_frames = 111679;
                num_cones = 1941;
                num_rgcs = 216;
                tc_length = 6;
                midget_window_size = 20;
                parasol_window_size = 45;
                save_path = ['~/Desktop/kiwi/'];                    
    otherwise
        disp('unknow data_name')
end

% load movie_xml_path & other info
datarun = load_index(datarun);

% load spikes times, trigger times, params, stas and cone info
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, datarun.names.nickname);
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

if p.Results.sta_summaries
    datarun = get_sta_summaries(datarun, {1,2,3,4}, 'keep_stas', false);
end

% load cone weights matrix
load([single_cone_path datarun.names.nickname '/Wc.mat'])
  
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


cone_info.Wc = Wc;
cone_info.cone_inputs = cone_inputs;
cone_info.spike_rate = spike_rate;
cone_info.spike_times = spike_times;
cone_info.rgc_tcs = rgc_tcs;
cone_info.rgc_ids = rgc_ids;
cone_info.cone_ids = cone_ids;


