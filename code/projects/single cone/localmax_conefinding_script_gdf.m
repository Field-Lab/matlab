%% load data
%datarun = load_data('2012-09-24-5/data001');

datarun = load_data('/Users/gdfield/Analysis/2011-10-25-5/data001/data001');
datarun = load_params(datarun);
info(datarun);

% set some piece parameters
datarun.piece.rig = 'A';
datarun.piece.optical_path_direction = 'below';
datarun.piece.display = 'oled1';
 
 
%% choose cells to work with
cellspec = {1,2,3,4,5};
datarun = load_sta(datarun, 'load_sta', []);
datarun = set_polarities(datarun);


%% Calculate STA summaries
marks_params.robust_std_method = 6;
marks_params.thresh = 4;
datarun = get_sta_summaries(datarun, cellspec, 'marks_params', marks_params);


%% Threshold each STA, scale by SNR, and add together

% get indices for cells of interest.
temp_indices = get_cell_indices(datarun, cellspec);
numRGCs = length(temp_indices);
 
% set threshold for significant stixels
stixel_threshold = 5;
 
% initialze a big matrix
normalized_RFs = zeros(numRGCs, prod(size(datarun.stas.rfs{temp_indices(1)})));
 
% process one cell at a time
for cc = 1:numRGCs
    % get cell ID and RF
    temp_cell_ID = datarun.cell_ids(temp_indices(cc));
    temp_rf = get_rf(datarun, temp_cell_ID);
    
    % normalize by the SNR
    temp_N = robust_std(reshape(temp_rf, 1,[]), 5);     % estimate the noise
    sig_stix = significant_stixels(temp_rf, 'thresh', stixel_threshold);
    temp_S = mean(temp_rf(sig_stix));     % estimate the signal
    temp_SNR = temp_S/temp_N; % estimate the SNR
    
    % threshold the RF and scale by SNR
    thresh_rf = zeros(size(temp_rf));
    thresh_rf(sig_stix) = temp_rf(sig_stix) * temp_SNR;
    normalized_RFs(cc,:) = reshape(thresh_rf, 1,[]);
    
    % Save for refinement stage
    
    all_sig_stix(:,cc) = sig_stix(:);
end

summed_rfs = reshape(sum(normalized_RFs), size(temp_rf));


%% Get the local maxima on the above image
temp_cone_locations = find_local_maxima(summed_rfs, 'radius', 1, 'thresh', 0, 'return', 'coordinates');

% plot cones on sensitivity image above
figure
imagesc(summed_rfs); colormap gray; hold on
plot(temp_cone_locations(:,1), temp_cone_locations(:,2), 'yo');
hold off
    
% For checking
cone_locations = temp_cone_locations;


%% Check how well we did with a given cell type, multiple figures
datarun = get_sta_fits_from_vision(datarun);

for cellid = datarun.cell_types{4}.cell_ids(1:19)
    figure();
    
    sanesubplot(2, 1, [1 1]);
    plot_rf(datarun, cellid);
    hold on;
    plot(cone_locations(:,1), cone_locations(:,2), 'yo');
    autozoom_to_fit(datarun, cellid, 5, 1, 1);
 
    sanesubplot(2, 1, [2 1]);
    imagesc(summed_rfs); colormap gray; axis equal tight ij
    hold on;
    plot(cone_locations(:,1), cone_locations(:,2), 'yo');
    autozoom_to_fit(datarun, cellid, 5, 1, 1);
end


%% Check how well we did with a given cell type, single panel
datarun = get_sta_fits_from_vision(datarun);
cellids = datarun.cell_types{4}.cell_ids(19:27);

rows = floor(sqrt(length(cellids)));
cols = ceil(length(cellids) / rows);
figure();
for i = 1:length(cellids)
    cellid = cellids(i);
    
    subplot(rows, cols, i);
    plot_rf(datarun, cellid, 'title', false);
    hold on;
    plot(cone_locations(:,1), cone_locations(:,2), 'yo');
    autozoom_to_fit(datarun, cellid, 5, 1, 1);
    title(num2str(cellid));
end


%% Now just need to define rough cone regions so we know which RGC's RFs contribute to each cone fit
mosaic = make_mosaic_struct(temp_cone_locations);

[V,C] = voronoin(temp_cone_locations);
% scale = datarun.stimulus.stixel_width;
masks = make_voronoi_masks(V, C, datarun.stimulus.field_width, datarun.stimulus.field_height);
clear V C scale;

near_cones_list = delaunay_neighbors(mosaic.delaunay);
neighbor_min_dist = 0;
self_max_dist = median(mosaic.nnd);
spaced = space_voronoi_masks(masks, temp_cone_locations, near_cones_list, neighbor_min_dist, self_max_dist);
clear near_cones_list neighbor_min_dist self_max_dist;

cones_labeled = index_masks(spaced, num2cell(1:length(spaced)));


%%
classify_params.algorithm = 'k means, EM';

if strcmpi(datarun.stimulus.independent, 'nil')
    % expand RFs to occupy all color channels
    for cc =1:length(datarun.cell_ids)
        if ~isempty(datarun.stas.rfs{cc})
            temp = datarun.stas.rfs{cc};
            datarun.stas.rfs{cc}(:,:,2) = temp;
            datarun.stas.rfs{cc}(:,:,3) = temp;
        end
    end
    clear temp
    
    % chose a phony-baloney cone classification algorithm
    classify_params.algorithm = 'nearest line';
end
clear cc


%% Refine cone center locations

cellids = datarun.cell_ids(temp_indices);

% Not working yet; need to get cones_labeled; used to be by running bwlabel
% on the sensitivity surface, but this probably needs to match the cones
% that are found as "initial centers" so there is more work to bring these
% into line.

%cone_rf_params.centers = 'fit'; % fit a gaussian
cone_rf_params.centers = 'com'; % find the center of mass
%cone_rf_params.centers = 'fixed';  % use the initial conditions
%cone_rf_params.centers = 'fit neighborhood';                % fit all cones in a small neighborhood
% cone_rf_params.sensitivity.strength = marks_params.strength; % how to combine RGB values of the RF
% cone_rf_params.sensitivity.filter = [];                     % how to filter the RF before looking for significant stixels
cone_rf_params.verbose = true;

if strcmpi(datarun.stimulus.independent, 't')
    cone_rf_params.single_cone_figure  = 105;
    cone_rf_params.cone_weights_figure = 106;
    
    cone_rf_params.cone_remap.fcn = @(x)([x(:,1)./x(:,2) x(:,3)./x(:,2)]);
    cone_rf_params.cone_remap.x_caption = 'red/green';
    cone_rf_params.cone_remap.y_caption = 'blue/green';
    
    cone_rf_params.regress = 1;
    cone_rf_params.combine = 'sum';
end


% FIXME: This presumably needs to be updated for OLED pixel sizes
% Cone kernel parameters; this shape is fitted to each cone
cone_kernel_params.type = 'dog';
cone_kernel_params.center_radius = 0.75;
cone_kernel_params.surround_scale = 0;
cone_kernel_params.surround_radius = 1;

[cone_rgb, cone_spatial_profiles, cone_ideal_shapes, cone_rfs, cone_ids] =                   ...
                        summarize_cone_rfs(datarun, cellids, temp_cone_locations,   ...
                        cones_labeled, all_sig_stix, cone_kernel_params, cone_rf_params);                    

                    
cone_centers = zeros(length(cone_ideal_shapes),2);
cone_found = true(length(cone_ideal_shapes),1);
for cc = 1:length(cone_ideal_shapes)
    if isempty(cone_ideal_shapes{cc}), 
        cone_centers(cc,1:2) = temp_cone_locations(cc,1:2);
        cone_found(cc) = false;
        continue; 
    end
    
    cone_centers(cc,1:2) = cone_ideal_shapes{cc}.center;
end

% For checking
cone_locations = cone_centers;


%% Check which cones were not found in refinement overall
figure
imagesc(summed_rfs); colormap gray; hold on
plot(cone_centers(cone_found,1),  cone_centers(cone_found,2),  'yo');
plot(cone_centers(~cone_found,1), cone_centers(~cone_found,2), 'ro');



%% Check how well we did with a given cell type, single panel, mark cones that were not found in refinement
datarun = get_sta_fits_from_vision(datarun);
cellids = datarun.cell_types{4}.cell_ids(13:24);

rows = floor(sqrt(length(cellids)));
cols = ceil(length(cellids) / rows);
figure();
for i = 1:length(cellids)
    cellid = cellids(i);
    
    subplot(rows, cols, i);
    plot_rf(datarun, cellid, 'title', false);
    hold on;    
    plot(cone_centers(cone_found,1),  cone_centers(cone_found,2),  'yo');
    plot(cone_centers(~cone_found,1), cone_centers(~cone_found,2), 'ro');
    autozoom_to_fit(datarun, cellid, 5, 1, 1);
    title(num2str(cellid));
end


%% Cut out cones that weren't found in refinement
% If things are working correctly, this should turn out to be only cones
% around the edges, where probably the issue is that the voronoi regions
% are incomplete and thus the masks are empty/truncated.  Could be fixed in
% part by fixing the Voronoi masks, but we don't care much about those edge
% cones usually.  For now punt and just cut them.

ncones = sum(cone_found);
cone_centers          = cone_centers(cone_found,:);
cone_rgb              = cone_rgb(cone_found,:);
cone_spatial_profiles = cone_spatial_profiles(:,cone_found);
cone_ideal_shapes     = cone_ideal_shapes(cone_found);
cone_rfs              = cone_rfs(cone_found);
cone_ids = 1:ncones;


%% Color classify cones
if strcmpi(datarun.stimulus.independent, 't')
    % classify cones
    [cone_types,likelihood_ratios,classify_extras] = ...
        classify_cones(cone_rgb, cone_rgb_expected(datarun),classify_params);
else
    cone_types = repmat('U', ncones, 1);
    classify_extras.cone_rgb_observed = struct('U', [1 1 1]);
end


%% Identify cone sampling by each RGC

% free up memory prior to this computation
%clear all_sig_stixels cone_ideal_shapes

[cone_weights, regr_cell_ids,Wc]= extract_cone_weights(datarun, cellspec,...
    cone_spatial_profiles, cone_types, classify_extras.cone_rgb_observed);


%% Fit center and surround gaussians to each RGC's cone RF
clear cone_info
cone_info.cone_weights = cone_weights;
cone_info.cone_centers = cone_centers;
cone_info.cell_ids = regr_cell_ids;

rf_cone_fits = fit_cone_rfs(datarun, cellspec, 'cone_info', cone_info, 'fit_radius', 150, 'verbose', 1, 'foa_2d', []);


%% Save stuff out
addname = ['localmaxgdf-st' num2str(stixel_threshold)];
savepath = sprintf('%s%s-%s/', single_cone_path, datarun.names.short_name, addname);
mkdir(savepath)

all_data.cone_weights = cone_weights;
all_data.cone_types = cone_types;
all_data.cone_centers = cone_centers;
all_data.rf_cone_fits = rf_cone_fits;

% Need to have even if it's bullshit
all_data.cone_rgb = cone_rgb; 
if isfield(classify_extras,'types_kmeans')
    all_data.cone_types_kmeans = classify_extras.types_kmeans;
    all_data.cone_types_em = classify_extras.types_em;
    all_data.cone_likelihoods = likelihood_ratios;
else
    all_data.cone_types_kmeans = cone_types;
    all_data.cone_types_em = cone_types;
    all_data.cone_likelihoods = ones(length(cone_types),1);
end

export_single_cone_data(datarun, cellspec, all_data, savepath);
datarun = load_cones(datarun, addname);

save([savepath 'Wc'],'Wc');

%FIXME: A few more data structs could be saved out, the parameters used is
%traditionally saved out, some figs could be saved out.


%% Save out stuff for Jeremy Freeman's analysis
datarun = load_neurons(datarun);
datarun = conepreprocess_wnm(datarun, 'cone_data_ind', addname);
conepreprocess_save(datarun, 'cone_data_ind', addname);