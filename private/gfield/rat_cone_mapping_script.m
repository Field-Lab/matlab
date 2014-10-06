
datarun = load_data('2011-09-23-0', 'rf-4');

datarun = load_index(datarun);
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun.stimulus = guess_stimulus(datarun.stimulus);
datarun = get_sta_summaries(datarun,'all','verbose',1,'keep_stas',0,'keep_rfs',1,'fig_or_axes',1,...
    'marks_params',struct('strength','vector length','filter',fspecial('gauss',15,0.7),'thresh',5));


datarun.stimulus.x_start = 0;
datarun.stimulus.x_end = 200;

datarun.stimulus.y_start = 0;
datarun.stimulus.y_end = 200;

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


start_time = 0;
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
bcf_params.cone_finding_cell_spec = {8,9};

% which cell types to identify the sampling of
bcf_params.regression_cell_spec = {8,9};
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
bcf_params.kernel_spacing = 1/6; bcf_params.kernel_radii = 0.5*[1 1 1];

% distance prior
bcf_params.C_C = [2.3 2.6]; 

% central square
% whole thing
bcf_params.relevant_region_x = [1 200-10]; bcf_params.relevant_region_y = [1 200-10];



%%
bcf = bayesian_cone_finding_loop(datarun,bcf_params);


% make some false cone loations and types if needed
datarun.cones.centers = [];
datarun.cones.types = [];

choose_magic_number(datarun,bcf,bcf_params);

save_bayesian_cones(datarun,bcf,bcf_params, 10, '10-0001-rat',false,[]);
