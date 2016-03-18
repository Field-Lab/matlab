%% Input here

% main parameters
datarun = load_data('/Volumes/Acquisition/Analysis/2016-03-17-2/data001/data001');
movie_descr = 'BW-3-6-0.48-11111-265x200-60.35.xml';
cell_types = {1,2,3,4,5}; 
% load('/Volumes/Analysis/2016-02-17-4/data001_denoised_sta', 'denoised_sta'); % if present

% additionbal parameters
datarun.names.nickname = '';
datarun.piece.rig = 'A'; % 'A' or 'B'
datarun.piece.optical_path_direction = 'below'; %'below' or 'above'
datarun.piece.display = 'oled1'; % 'oled' or 'crt' with 1 for A, 2 for B


%% load stuff
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

if exist('denoised_sta', 'var')
    datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
    for i=1:length(datarun.cell_ids)
        tmp = denoised_sta(:,:,:,i);
        datarun.stas.stas{i} = permute(tmp, [1 2 4 3]);
    end
else
    datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
end

for i=1:length(datarun.cell_ids)
    tmp = datarun.stas.stas{i}-repmat(denoised_sta,1,1,1,6);
    datarun.stas.stas{i} = tmp;
end

datarun = set_polarities(datarun);
% find movie, check if exists
movie_spec = fullfile('/Volumes/Analysis/stimuli/white-noise-xml/', movie_descr);
independent = strcmpi(datarun.stimulus.independent, 't');
field_width = datarun.stimulus.field_width;
field_height = datarun.stimulus.field_height;

cell_inds = get_cell_indices(datarun, cell_types);

%% get sta info - summary and static nonlinearities

% robust_std_method is 1 to match old implementation.  Set to 3,5 for some speedup.
disp('Loading sta summaries and static NLs...')
datarun = get_sta_summaries(datarun, cell_types, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
        'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
        'thresh',5,'robust_std_method',1));
    
% calculate static nonlinearities
datarun = load_java_movie(datarun, movie_spec);
datarun = get_snls(datarun, datarun.cell_ids(cell_inds),...
    'frames',-2:0,'start_time',0,'stimuli',3000,'new',true); % most memory and time-consuming part

%% initial cone estimation, prior set up (local max)
stixel_threshold = 4.5;
thresh = 15;
radius = 1;

combine_stixels = 'sum';
% combine_stixels = 'max';
[surface,all_sig_stixels,spatial_cell_ids] = compute_spatial_sensitivity(datarun, ...
    cell_types,'verbose', true, 'combine_stixels', ...
    combine_stixels, 'selection_params', struct('type', 'max', 'thresh', stixel_threshold));
surface(surface < thresh) = 0;

local_maxima = find_cones_in_rf(surface, 'filter', [], ...
    'selection', struct('type', 'max', 'thresh', 0, 'radius', radius));

cones_labeled = bwlabel(local_maxima,8);
ncones = max(cones_labeled(:));
cone_centers = zeros(ncones,2);
for nn = 1:ncones
    [cone_centers(nn,1),cone_centers(nn,2)] = ait_centroid(cones_labeled == nn);
end
clear nn

min_dist = pdist2(cone_centers,cone_centers);
min_dist(min_dist==0) = 100;
min_dist = min(min_dist);

figure
subplot(1,2,1)
imagesc(surface.^0.5)
colormap gray
hold on; plot(cone_centers(:,1), cone_centers(:,2), 'yo')
subplot(1,2,2)
hist(min_dist,0:0.5:8)
axis([0 8 0 Inf])
title('Nearest Neighbor distance for cones. SET PRIOR!')

myPrior=[2.75,3.25];

%% set bcf params

bcf_params = struct;
bcf_params.kernel_plot_colors = ('g')';
bcf_params.C_C = myPrior; % distance prior
% bcf_params.kernel_colors = cone_rgb_expected(datarun);

% figure
bcf_params.plot_fig = 20;
bcf_params.cones_fig = 21;
bcf_params.dll_fig = 22;

% cell spec
bcf_params.cone_finding_cell_spec = cell_types;

% size of subsets
bcf_params.padding_x = 5; % was 5 before
bcf_params.padding_y = 5; % was 5 before
% bcf_params.roi_x_size = 2*bcf_params.padding_x + 10; % was +10 before. in general, should be +(pad size)*2
% bcf_params.roi_y_size = 2*bcf_params.padding_y + 10; % was +10 before

bcf_params.roi_x_size = 20; % was +10 before. in general, should be +(pad size)*2
bcf_params.roi_y_size = 20; % was +10 before

% don't recompute if not necessary
bcf_params.new_W = 0;
bcf_params.new_STAs = 0;

% iterations per patch
bcf_params.num_iter = 100;

% radius of relevance
% sets which stixels (around the marks) are considered relevant in each STA
bcf_params.rel_radius = 4;

% start time of the datarun
% bcf_params.start_time = tic;

% arbitrary scale factor
bcf_params.magic_number = 1;

% cone density prior
bcf_params.q = 0.05;

% kernel spacing (in pixels) and radius
bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

bcf_params.relevant_region_x = [1 field_width-2*bcf_params.padding_x]; bcf_params.relevant_region_y = [1 field_height-2*bcf_params.padding_y];  % x == width, y == height??

%%
% note start
start_time_all = clock;
bcf = bcf_loop(datarun,bcf_params);
fprintf('\nfinished loop in %0.1f sec\n\n\n',etime(clock,start_time_all))


%%
datarun.cones.centers = [];
datarun.cones.types = [];
choose_magic_number(datarun,bcf,bcf_params);


%%
magic_number = 20;
path2save = ['/Volumes/Analysis/2016-03-17-2/cone_data/data001/denoised_bayes-msf_', int2str(magic_number),'/'];
mkdir(path2save)
keep_indices = (bcf.all_added_cones(:,7) + magic_number*bcf.all_added_cones(:,8)) > 0;
cone_centers = bcf.all_added_cones(keep_indices,[4 5]);
datarun.cones.centers = cone_centers;
save([path2save, 'cones'], 'cone_centers', 'bcf', 'bcf_params');

%%
% load('/Volumes/Analysis/2016-02-17-4/cone_data/data001/denoised_bayes-msf_70/cones.mat', 'cone_centers')
% datarun.cones.centers = cone_centers;
datarun = make_mosaic_struct(datarun);
datarun = make_voronoi_masks(datarun);
masks = datarun.cones.mosaic.voronoi_masks;

min_neighbor_dist = 2;
max_self_dist = 4;
spaced = space_voronoi_masks(datarun, min_neighbor_dist, max_self_dist);
cone_map = index_masks(spaced, num2cell(1:length(masks)));

figure
imagesc(cone_map)

dlmwrite(['/Volumes/Analysis/2016-03-17-2/cone_data/data001/denoised_bayes-msf_', int2str(magic_number),'/map_data001.txt'], cone_map, 'delimiter', '\t', 'newline', 'pc')

a = dlmread(['/Volumes/Analysis/2016-03-17-2/cone_data/data001/denoised_bayes-msf_', int2str(magic_number),'/map_data001.txt']);
figure;
imagesc(a)


%%

tmp = load('/Volumes/Analysis/2016-02-17-4/cone_data/data001/denoised_bayes-msf_70/cones.mat', 'cone_centers')
tmp1 = load('/Volumes/Analysis/2016-02-17-4/cone_data/data001/denoised_bayes-msf_120/cones.mat', 'cone_centers')

figure
plot(tmp.cone_centers(:,1),tmp.cone_centers(:,2), 'xr' )
hold on
plot(tmp1.cone_centers(:,1),tmp1.cone_centers(:,2), '+g' )

