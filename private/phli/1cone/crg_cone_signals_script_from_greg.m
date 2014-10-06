%% LOAD DATA 
% get datarun with cone map
opt=struct('verbose',1,'load_params',1,'load_neurons',1);
datarun{1} = load_data('2010-03-05-2/data013', opt);
datarun{1} = load_sta(datarun{1}, 'load_sta', 'all', 'save_sta', false);
datarun{1} = get_sta_summaries(datarun{1}, {1,2,3,4,5}, 'keep_stas', false, 'keep_rfs', false, 'marks_params', struct('robust_std_method', 6));
%datarun{1} = import_single_cone_data(datarun{1}, 'path');

% get grating datarun
datarun{2} = load_data('2010-03-05-2/data016-from-data013', opt);

% map cells
datarun=map_cell_types(datarun);

% parse stimulus
datarun{2}.names.stimulus_path = '/snle/analysis/2010-03-05-2/stimuli/s16';
datarun{2}=load_stim(datarun{2}); 

%% Cone information


% load cone info
%if ~isfield(datarun, 'cones') || isempty(datarun.cones.centers)
%    datarun = load_cones(datarun, opts.cone_data_ind);
%end
conepath = [single_cone_path '/apple/'];
load(fullfile(conepath, 'Wc.mat'));

field_width = datarun{1}.stimulus.field_width;
field_height = datarun{1}.stimulus.field_height;

datarun{1} = import_single_cone_data(datarun{1},conepath);

%% Grating stimuli
% set invariant stimulus parameters
temp_spec.x_start = 0;
temp_spec.y_start = 0;
temp_spec.x_end = 640;
temp_spec.y_end = 480;
temp_spec.orientation = 0;
grating_duration = 8;
frame_rate = 120;

for stim_combo = 1:length(datarun{2}.stimulus.combinations)

    temp_spec.spatial_period = datarun{2}.stimulus.combinations(stim_combo).SPATIAL_PERIOD;
    temp_spec.temporal_period = datarun{2}.stimulus.combinations(stim_combo).TEMPORAL_PERIOD;
    temp_spec.spatial_phase = datarun{2}.stimulus.combinations(stim_combo).SPATIAL_PHASE;


    [temp_frame, temp_tscale] = calc_reversing_grating_frame_intensities(temp_spec);
    
    stimulus_set(stim_combo).max_frame = temp_frame;
    stimulus_set(stim_combo).tscale = repmat(temp_tscale, 1, grating_duration * frame_rate / temp_spec.temporal_period);
    
end

%% compute cone generator sigs in response to stimulus

cone_inputs = zeros(size(Wc,2), length(stimulus_set(1).tscale), length(stimulus_set));
% align cones with stimulus and compute amplitude modulation in time.
for stim_combo = 1:length(stimulus_set)
    max_frame = reshape(repmat(stimulus_set(stim_combo).max_frame(81:400,161:480), [1,1,3]),1,[]);
    stim_space_by_time = max_frame' * stimulus_set(stim_combo).tscale;
    cone_inputs(:,:,stim_combo) = Wc' * stim_space_by_time;
end

datarun{2}.cone_inputs = cone_inputs;

%%
figure
load ~/Desktop/gratingCode/dat-826
g = getGratingResp(datarun{2},dat);
plotGratings(g,datarun{2},dat,'cycle');

%%
print(3, '/Users/peterli/Desktop/off-par-826.pdf', '-dpdf')

%dat.rgcId -- the cell id for the cell you want to look at
%dat.coneIds -- the cone numbers for this cell (these are just indices into the full list of cones, i.e. the first dimension of the cone_input matrix)
%dat.locs_c -- the xy locations for those cones in the piece (for plotting)
%dat.xRange -- range of x values (for setting the plot axes)
%dat.yRange -- range of y values (for setting the plot axes)

    
    
    
    








