starttime = tic;

datarun = load_data(fullfile(server_path(), '2012-09-13-2/data009/data009'));

datarun.piece.rig = 'A';
datarun.piece.optical_path_direction = 'below';
datarun.piece.display = 'oled1';


%% Load cell types, check a few other spots before resorting to params
datarun = load_txt_cell_types(datarun, 'localmaxconefinding');
datarun = load_txt_cell_types(datarun, 'conefinding');
datarun = load_txt_cell_types(datarun, '1coneclassification');
datarun = load_params(datarun);
info(datarun);


%%
cellspec = {1,2,3,4,5};


%%
datarun = load_sta(datarun, 'load_sta', []);
datarun = set_polarities(datarun);


%%
marks_params.robust_std_method = 6;
marks_params.thresh = 4;
if strcmpi(datarun.stimulus.independent, 't')
    rgb_expected = cone_rgb_expected(datarun);
    marks_params.strength = {'inner or', [rgb_expected.L;    ...
                             rgb_expected.M;    ...
                             rgb_expected.S]};
else
    marks_params.strength = 'vector length'; 
end

datarun = get_sta_summaries(datarun, cellspec, 'marks_params', marks_params);


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


%% Get spatial sensitivity from all RFs combined

% sig_stixel threshold in sigma
stixel_threshold = 4.5;

disp('Compute spatial sensitivity...');
combine_stixels = 'sum';
% combine_stixels = 'max';
[spatial_sensitivity,all_sig_stixels,spatial_cell_ids] = compute_spatial_sensitivity(datarun, cellspec, 'verbose', true, 'combine_stixels', combine_stixels, 'selection_params', struct('type', 'max', 'thresh', stixel_threshold));
hist(spatial_sensitivity(:), 100)


%% Simple thresh cutoff
% Statistical thresh calculation below is also a simple threshold, but the cutoff value is calculated from the sigma requested
thresh = 25;
thresh_spatial_sensitivity = spatial_sensitivity;
thresh_spatial_sensitivity(thresh_spatial_sensitivity < thresh) = 0;


%% Get initial estimate of cone centers
surface = thresh_spatial_sensitivity;
radius = 1;

disp('Initial cone estimates...');
tic
local_maxima = find_cones_in_rf(surface, 'filter', [], 'selection', struct('type', 'max', 'thresh', 0, 'radius', radius));
cones_labeled = bwlabel(local_maxima,8);
ncones = max(cones_labeled(:));

initial_cone_centers = zeros(ncones,2);
for nn = 1:ncones
    [initial_cone_centers(nn,1),initial_cone_centers(nn,2)] = ait_centroid(cones_labeled == nn);
end
clear nn
toc

% If we're planning to skip the refinement step for now...
cone_centers = initial_cone_centers;

imagesc(surface.^0.5)
colormap gray
hold on; plot(initial_cone_centers(:,1), initial_cone_centers(:,2), 'yo')


%% Check how well we did with a given cell type, single panel
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
    plot(cone_centers(:,1), cone_centers(:,2), 'yo');
    autozoom_to_fit(datarun, cellid, 5, 1, 1);
    title(num2str(cellid));
end


%% Refine cone center locations

%cone_rf_params.centers = 'fit'; % fit a gaussian
cone_rf_params.centers = 'com'; % find the center of mass
%cone_rf_params.centers = 'fixed';  % use the initial conditions
%cone_rf_params.centers = 'fit neighborhood';                % fit all cones in a small neighborhood
cone_rf_params.sensitivity.strength = marks_params.strength; % how to combine RGB values of the RF
cone_rf_params.sensitivity.filter = [];                     % how to filter the RF before looking for significant stixels
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

[cone_rgb, cone_spatial_profiles, cone_ideal_shapes, cone_rfs, cone_ids] =                      ...
                        summarize_cone_rfs(datarun, spatial_cell_ids, initial_cone_centers,     ...
                        cones_labeled, all_sig_stixels, cone_kernel_params, cone_rf_params);                    

                    
cone_centers = zeros(length(cone_ideal_shapes),2);
for cc = 1:length(cone_ideal_shapes)
    cone_centers(cc,1:2) = cone_ideal_shapes{cc}.center;
end

% imagesc(spatial_sensitivity)
% colormap gray
% hold on; plot(initial_cone_centers(:,1), initial_cone_centers(:,2), 'yo')
% plot(cone_centers(:,1), cone_centers(:,2), 'bo')


%% Check how well we did with a given cell type
datarun = get_sta_fits_from_vision(datarun);

for cellid = datarun.cell_types{4}.cell_ids(1:20)
    figure();
    
    sanesubplot(2, 1, [1 1]);
    plot_rf(datarun, cellid);
    hold on;
    plot(cone_centers(:,1), cone_centers(:,2), 'yo');
%     plot(cone_centers(cone_types == 'L',1), cone_centers(cone_types == 'L',2), 'ro');
%     plot(cone_centers(cone_types == 'M',1), cone_centers(cone_types == 'M',2), 'go');
%     plot(cone_centers(cone_types == 'S',1), cone_centers(cone_types == 'S',2), 'bo');
    autozoom_to_fit(datarun, cellid, 5, 1, 1);

    sanesubplot(2, 1, [2 1]);
    imagesc(surface); colormap gray; axis equal tight ij
    hold on;
    plot(cone_centers(:,1), cone_centers(:,2), 'yo');
    autozoom_to_fit(datarun, cellid, 5, 1, 1);
    
%     sanesubplot(2, 2, [1 2]);
%     plot_rf(datarun, cellid);
%     hold on;
%     plot(datarun2.cones.centers(:,1), datarun2.cones.centers(:,2), 'yo');
%     autozoom_to_fit(datarun, cellid);
%     
%     sanesubplot(2, 2, [2 2]);
%     imagesc(bcf.dll.^(1/3)); colormap gray; axis equal tight ij
%     hold on;
%     plot(datarun2.cones.centers(:,1).*dll_scale, datarun2.cones.centers(:,2).*dll_scale, 'yo');
%     autozoom_to_fit(datarun, cellid, 5, dll_scale);
end

%% Check how well we did versus Bayesian
datarun = get_sta_fits_from_vision(datarun);

datarun2 = datarun;
datarun2 = load_cones(datarun2, 'bayes-msf_10.00');


% Load Bayesian DLL and scale
load([single_cone_path datarun2.names.cones '/results.mat']); % bcf
load([single_cone_path datarun2.names.cones '/parameters.mat']); % bcf_params
dll_scale = 1 ./ bcf_params.kernel_spacing;

for cellid = datarun.cell_types{4}.cell_ids(1:20)
    figure();
    
    sanesubplot(2, 2, [1 1]);
    plot_rf(datarun, cellid);
    hold on;
    plot(cone_centers(:,1), cone_centers(:,2), 'yo');
%     plot(cone_centers(cone_types == 'L',1), cone_centers(cone_types == 'L',2), 'ro');
%     plot(cone_centers(cone_types == 'M',1), cone_centers(cone_types == 'M',2), 'go');
%     plot(cone_centers(cone_types == 'S',1), cone_centers(cone_types == 'S',2), 'bo');
    autozoom_to_fit(datarun, cellid, 5, 1, 1);

    sanesubplot(2, 2, [2 1]);
    imagesc(surface); colormap gray; axis equal tight ij
    hold on;
    plot(cone_centers(:,1), cone_centers(:,2), 'yo');
    autozoom_to_fit(datarun, cellid, 5, 1, 1);
    
    sanesubplot(2, 2, [1 2]);
    plot_rf(datarun, cellid);
    hold on;
    plot(datarun2.cones.centers(:,1), datarun2.cones.centers(:,2), 'yo');
    autozoom_to_fit(datarun, cellid, 5, 1, 1);
    
    sanesubplot(2, 2, [2 2]);
    imagesc(bcf.dll.^(1/3)); colormap gray; axis equal tight ij
    hold on;
    plot(datarun2.cones.centers(:,1).*dll_scale, datarun2.cones.centers(:,2).*dll_scale, 'yo');
    autozoom_to_fit(datarun, cellid, 5, dll_scale, 1);
end


%% Estimate density, and only include cones in which local density is close to average density
% 
% % parameters
% num_bins = 20;
% bin_width = 1;
% 
% % get local density
% [drp,bin_centers,drp_extras] = density_recovery_profile(cone_centers, num_bins, bin_width);
% 
% % identify ROI
% cone_roi = identify_cone_mosaic_roi(cone_centers,drp_extras.density, drp_extras.eff_rad);
% 
% 
% imagesc(spatial_sensitivity)
% colormap gray
% hold on;
% plot(cone_centers(:,1), cone_centers(:,2), 'bo')
% plot(cone_centers(cone_roi,1), cone_centers(cone_roi,2), 'w.')


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
addname = ['ath25-localmax-t' num2str(thresh)];
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

save([savepath 'Wc'],'Wc')

%FIXME: A few more data structs could be saved out, the parameters used is
%traditionally saved out, some figs could be saved out.


%% Undo fake RGB RF
if strcmpi(datarun.stimulus.independent, 'nil')
    datarun.stas.rfs = cellfun(@(M)(M(:,:,1)), datarun.stas.rfs, 'UniformOutput', false);
end


%% Save out stuff for Jeremy Freeman's analysis
datarun = load_neurons(datarun);
datarun = conepreprocess_wnm(datarun);
conepreprocess_save(datarun, 'cone_data_ind', addname, 'cell_types', cellspec);

%%
toc(starttime);