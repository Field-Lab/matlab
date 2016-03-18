%% Input here

% main parameters
datarun = load_data('/Volumes/Acquisition/Analysis/2016-03-17-2/data001/data001');
movie_descr = 'BW-3-6-0.48-11111-265x200-60.35.xml';
cell_types = {1,2,3,4,5}; 
wid = 265;
hei = 200;
% load('/Volumes/Analysis/2016-02-17-4/data001_denoised_sta', 'denoised_sta'); % if present

% additionbal parameters
datarun.names.nickname = '';
datarun.piece.rig = 'A'; % 'A' or 'B'
datarun.piece.optical_path_direction = 'below'; %'below' or 'above'
datarun.piece.display = 'oled1'; % 'oled' or 'crt' with 1 for A, 2 for B


%% load stuff
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
datarun = set_polarities(datarun);
cell_inds = get_cell_indices(datarun, cell_types);

% denoise
% collect STAs
all_sta = zeros(hei,wid,length(datarun.cell_ids));
for i = 1:length(datarun.cell_ids)
    tmp = datarun.stas.stas{i};
    t = zeros(1,6);
    for j=1:6
        sta=squeeze(tmp(:,:,:,j));
        t(j) = max(abs(sta(:)));
    end
    [~, frame] = max(t);
    sta=tmp(:,:,:,frame);
    all_sta(:,:,i)=double(squeeze(sta));
end
noise_sta = mean(all_sta(:,:,cell_inds),3);
for i=1:length(datarun.cell_ids)
    tmp = datarun.stas.stas{i}-repmat(noise_sta,1,1,1,6);
    datarun.stas.stas{i} = tmp;
end

% find movie, check if exists
movie_spec = fullfile('/Volumes/Analysis/stimuli/white-noise-xml/', movie_descr);
independent = strcmpi(datarun.stimulus.independent, 't');
field_width = datarun.stimulus.field_width;
field_height = datarun.stimulus.field_height;


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
% 
% figure
% subplot(1,2,1)
% imagesc(surface.^0.5)
% colormap gray
% hold on; plot(cone_centers(:,1), cone_centers(:,2), 'yo')
% subplot(1,2,2)
% hist(min_dist,0:0.5:8)
% axis([0 8 0 Inf])
% title('Nearest Neighbor distance for cones. SET PRIOR!')

myPrior=[2.75,3.25];

%% set bcf params

bcf_params = struct;
bcf_params.kernel_plot_colors = ('g')';
bcf_params.C_C = myPrior; % distance prior

% cell spec
bcf_params.cone_finding_cell_spec = cell_types;

% size of subsets
bcf_params.padding_x = 5; % was 5 before
bcf_params.padding_y = 5; % was 5 before
% bcf_params.roi_x_size = 2*bcf_params.padding_x + 10; % was +10 before. in general, should be +(pad size)*2
% bcf_params.roi_y_size = 2*bcf_params.padding_y + 10; % was +10 before

bcf_params.roi_x_size = 15; % was +10 before. in general, should be +(pad size)*2
bcf_params.roi_y_size = 15; % was +10 before

% don't recompute if not necessary
bcf_params.new_W = 0;
bcf_params.new_STAs = 0;

% iterations per patch
bcf_params.num_iter = 100;

% radius of relevance
% sets which stixels (around the marks) are considered relevant in each STA
bcf_params.rel_radius = 4;

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

save('/Volumes/Analysis/2016-03-17-2/cone_data/data001/bcf.mat', 'bcf', 'bcf_params');
