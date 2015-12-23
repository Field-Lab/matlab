%% Input here
starttime = tic;
datarun = load_data(fullfile(server_path(), '2012-09-24-5/data001/data001'));

datarun.names.nickname = '';
datarun.piece.rig = 'A';
datarun.piece.optical_path_direction = 'below';
datarun.piece.display = 'oled1';
extra_dirname_info = 'BW-2-5';


%%
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
%datarun.stimulus = guess_stimulus(datarun.stimulus);

% BW or RGB stimulus?
independent = strcmpi(datarun.stimulus.independent, 't');
field_width = datarun.stimulus.field_width;
field_height = datarun.stimulus.field_height;

info(datarun);


%%
cell_types = {1,2,3,4,5};

tic;
% robust_std_method is 1 to match old implementation.  Set to 3,5 for some speedup.
datarun = get_sta_summaries(datarun, cell_types, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
        'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
        'thresh',5,'robust_std_method',3));
toc



%%

% load java movie
% %datarun = load_java_movie(datarun);
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


%% set bcf params
% initialize struct
bcf_params = struct;

switch independent
    case false
        % Fake a RGB STA from BW
        for cc = 1:length(datarun.stas.rfs)
            temp_rf = datarun.stas.rfs{cc};
            temp_rf = repmat(temp_rf, [1 1 3]);
            datarun.stas.rfs{cc} = temp_rf;
        end

        % BW use green
        bcf_params.kernel_plot_colors = ('g')';
        bcf_params.C_C = [2.3 2.5]; % distance prior

    case true % RGB stimulus
        bcf_params.kernel_plot_colors = ('rgb')';
        bcf_params.LM_MM = [2.3 2.5];
        bcf_params.LM_S = [1.8 2.0];
        bcf_params.S_S = [3.8 4];
end

% figure
bcf_params.plot_fig = 20;
bcf_params.cones_fig = 21;
bcf_params.dll_fig = 22;

% cell spec
bcf_params.cone_finding_cell_spec = cell_types;

% which cell types to identify the sampling of
bcf_params.regression_cell_spec = cell_types;

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

% central square
% whole thing
bcf_params.relevant_region_x = [1 field_width-10]; bcf_params.relevant_region_y = [1 field_height-10];  % x == width, y == height??



%%
tic;
bcf = bayesian_cone_finding_loop(datarun,bcf_params);
toc

toc(starttime);


%%
% make some false cone locations and types if needed
datarun.cones.centers = [];
datarun.cones.types = [];
    
choose_magic_number(datarun,bcf,bcf_params);


%%
magic_number = 5;
save_bayesian_cones(datarun, bcf, bcf_params, magic_number, extra_dirname_info, false, [], 'fit_foa', [], 'robust_std_method', 5);
toc


%% Undo fake RGB RF
if ~independent
    datarun.stas.rfs = cellfun(@(M)(M(:,:,1)), datarun.stas.rfs, 'UniformOutput', false);
end


%% Save out stuff for Jeremy Freeman's analysis
datarun = load_neurons(datarun);
datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
conepreprocess_save(datarun, 'cone_data_ind', 'bayes');