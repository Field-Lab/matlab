
%datarun = load_data('/snle/acquisition/2010-09-24-1/data002/data002');
%datarun = load_data('/snle/acquisition/2010-09-24-1/data006/data006');
%datarun = load_data('/snle/acquisition/2010-09-24-1/data035/data035');

% path during experiment
datarun = load_data('/snle/acquisition/2010-09-24-1/data006/data006');

% path for testing
datarun = load_data('/snle/lab/Experiments/Array/Analysis/2010-09-24-1/streamed/data006/data006');


%datarun = load_index(datarun);
%datarun = load_data(datarun_spec);
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun.stimulus = guess_stimulus(datarun.stimulus);
datarun = get_sta_summaries(datarun,{1,2,3,4,5},'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct('strength','vector length','filter',fspecial('gauss',15,0.7),'thresh',5));


datarun.names.nickname = 'nickname';

% load java movie
%datarun = load_java_movie(datarun);
% should use thousands of stimuli, or all
start_time = 0;
%datarun = get_snls(datarun, {1,2,3,4,5},'frames',:,'start_time',start_time,'stimuli',10000,'new',false);
%save([single_cone_path 'saved/' datarun.names.nickname],'datarun','-v7.3')

%% put in fake SNL information
cell_indices = get_cell_indices(datarun, 'all');
fit_params.a = 1;
fit_params.b = 1;
fit_params.type = 'exp';
for cc = 1:length(cell_indices)
    datarun.stas.snls{cell_indices(cc)}.fit_params = fit_params;
end


%% If needed: fake a RGB STA from BW

for cc = 1:length(datarun.stas.rfs)
    temp_rf = datarun.stas.rfs{cc};
    temp_rf = repmat(temp_rf, [1 1 3]);
    datarun.stas.rfs{cc} = temp_rf;
end

%% set bcf params
% initialize struct
bcf_params = struct;

% MIGHT HAVE TO CHANGE KERNEL_PLOT_COLORS IF STIMULUS IS RGB
% kernel plot colors
bcf_params.kernel_plot_colors = ('g')';
%bcf_params.kernel_plot_colors = ('rgb')';

% figure
bcf_params.plot_fig = 20;
bcf_params.cones_fig = 21;
bcf_params.dll_fig = 22;

% cell spec
bcf_params.cone_finding_cell_spec = {1,2,3,4,5};

% which cell types to identify the sampling of
bcf_params.regression_cell_spec = {1,2,3,4,5};
%bcf_params.regression_cell_spec = 'all';

% size of subsets
bcf_params.padding_x = 5;
bcf_params.padding_y = 5;
bcf_params.roi_x_size = 2*bcf_params.padding_x + 10;
bcf_params.roi_y_size = 2*bcf_params.padding_y + 10;

% don't recompute if not necessary
bcf_params.new_W = 0;
bcf_params.new_STAs = 0;

% iterations per patch
bcf_params.num_iter = 100 ;

% radius of relevance
% sets which stixels (around the marks) are considered relevant in each STA
bcf_params.rel_radius = 4;

% kernel RGB
bcf_params.kernel_colors = cone_rgb_expected(datarun);


% start time of the datarun
bcf_params.start_time = start_time;

% arbitrary scale factor
bcf_params.magic_number = 1;

% cone density prior
bcf_params.q = 0.05;

% kernel spacing (in pixels) and radius
bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

% distance prior
bcf_params.C_C = [2.3 2.5]; 

% central square
% whole thing
bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];


%%

% Save the setup to disk to be loaded by parallel machines
%save ...

% Load the setup on each parallel machine
%load ...


%%


% Setup machine numbers
machines = struct();
machines.smokestack = 0;
machines.node67     = 1;

[~, this_machine] = system('hostname -s');
rois_index = machines.(this_machine);


bcfs{rois_index} = bayesian_cone_finding_split_loop(datarun, bcf_params, 'thread', rois_index, 'threads', 2);


% Save bcfs to disk
% save ...

%%

% Load all bcfs
% for ... load ...

bcf = combine_bcfs(bcfs);

%%

% make some false cone loations and types if needed
datarun.cones.centers = [];
datarun.cones.types = [];

choose_magic_number(datarun,bcf,bcf_params);

save_bayesian_cones(datarun,bcf,bcf_params, 15, '15-foo',false,[]);
