%% LOCAL MAX CONE FINDING
% This code is an attempt to simpify and streamline earlier versions of
% local-max cone finding.  
% :GDF 2013-02

%% load data

% CRT data

datarun = load_data('/Volumes/Analysis/2010-03-05-2/data001/data001');

% set some piece parameters
datarun.piece.rig = 'A';
datarun.piece.optical_path_direction = 'below';
datarun.piece.display = 'crt1';

% OLED data
datarun = load_data('/Volumes/Analysis/2012-09-13-2/data005/data005');
datarun =  load_data('/Volumes/Analysis/2011-10-25-5/data001/data001');
datarun = load_data(fullfile(server_path(), '2012-09-13-2/data009/data009'));

% set some piece parameters
datarun.piece.rig = 'A';
datarun.piece.optical_path_direction = 'below';
datarun.piece.display = 'oled1';

%%


datarun = load_params(datarun);
info(datarun);

gen_params.marks_thresh = 5; % sigmas: used to calculate the marks for RFs and TCs
gen_params.stixel_threshold = 5; % sigmas: used to calc marks for cell-wide sensitivity map
gen_params.cone_location_optimization_roi = 1; % ROI radius for calculating the COM for each cone
gen_params.cone_size_roi = 3; % ROI radius for estimating the mean cone size
gen_params.cone_upsampling = 10; % factor by which to upsample to cone RFs for interpolation when estimating their size
gen_params.cone_rgb_roi = 1; % ROI radius for estimating the rgb values assoicated with each cone.


%% choose cells to work with
cellspec = {1,2,3,4,5};
datarun = load_sta(datarun, 'load_sta', []);
datarun = set_polarities(datarun);


%% Calculate STA summaries
marks_params.robust_std_method = 6;
marks_params.thresh = gen_params.marks_thresh;
datarun = get_sta_summaries(datarun, cellspec, 'marks_params', marks_params);

% recompute for SBCs and their SNR tends to be lower
if strcmp(datarun.stimulus.independent, 't');
    sbc_marks_params.thresh = 4;
    sbc_marks_params.strength = {'inner',[0 0 1]};
    datarun = get_sta_summaries(datarun, {5}, 'marks_params', sbc_marks_params);
end

%% Threshold each STA, scale by SNR, and add together

% set threshold for significant stixels
stixel_threshold = gen_params.stixel_threshold; 

% get indices for cells of interest.
temp_indices = get_cell_indices(datarun, cellspec);
numRGCs = length(temp_indices);
 
% initialze a big matrix
normalized_RFs = zeros(numRGCs, prod([size(datarun.stas.rfs{temp_indices(1)},1), size(datarun.stas.rfs{temp_indices(1)},2)]));
all_sig_stix = zeros(prod([size(datarun.stas.rfs{temp_indices(1)},1), size(datarun.stas.rfs{temp_indices(1)},2)]), numRGCs);
% process one cell at a time
for cc = 1:numRGCs
    % get cell ID and RF
    temp_cell_ID = datarun.cell_ids(temp_indices(cc));
    temp_rf = get_rf(datarun, temp_cell_ID);        
    
    % check to see if stim is RGB
    if size(temp_rf, 3) == 3
       
        % check to see if SBC
        if ~isempty(find(datarun.cell_types{5}.cell_ids == temp_cell_ID))
            
            % get just the blue gun
            temp_rf = squeeze(temp_rf(:,:,3));
            
            % get sig stixel
            sig_stix = significant_stixels(temp_rf, 'thresh', 5);
            
%             figure(1)
%             imagesc(temp_rf); colormap gray; 
%             pause
%             imagesc(sig_stix)
%             pause
        else
            % otherwise sum all guns
            temp_rf = sum(temp_rf, 3);

            sig_stix = significant_stixels(temp_rf, 'thresh', stixel_threshold);
        end
    else
        % if already BW, calculate significant stixels
        sig_stix = significant_stixels(temp_rf, 'thresh', stixel_threshold);
    end
    
    if ~isempty(temp_rf)
        % normalize by the SNR
        temp_N = robust_std(reshape(temp_rf, 1,[]), 5);     % estimate the noise
        temp_S = mean(temp_rf(sig_stix));     % estimate the signal
        temp_SNR = temp_S/temp_N^2; % estimate the SNR
        
        % threshold the RF and scale by SNR
        thresh_rf = zeros(size(temp_rf));
        thresh_rf(sig_stix) = temp_rf(sig_stix) * temp_SNR;
        normalized_RFs(cc,:) = reshape(thresh_rf, 1,[]);
        
        % Save for refinement stage
        
        all_sig_stix(:,cc) = sig_stix(:);
    end
end

summed_rfs = reshape(sum(normalized_RFs), size(temp_rf));


%% Get the local maxima on the above image
% 
% load('/Users/alexth/Desktop/dll.mat')
% a=a(1:9:end,1:9:end,1);
% figure
% imagesc(a)
% figure
% imagesc(summed_rfs)
% cone_peaks = find_local_maxima(a, 'radius', 1, 'thresh', 0, 'return', 'indices');
% 
% figure
% imagesc(a); colormap gray; hold on
% plot(cone_peaks(:,2), cone_peaks(:,1), 'yo');
% hold off

cone_peaks = find_local_maxima(summed_rfs, 'radius', 1, 'thresh', 0, 'return', 'indices');

% plot cones on sensitivity image above
figure
imagesc(summed_rfs); colormap gray; hold on
plot(cone_peaks(:,2), cone_peaks(:,1), 'yo');
hold off


%% ath distance measure

for i=1:length(cone_peaks)
    
    D=pdist2(cone_peaks(i,:),cone_peaks);
    m(i)=min(D(D>0));
end

figure
hist(m,0:0.5:10)



%% Check how well we did with a given cell type, multiple figures
        datarun = get_sta_fits_from_vision(datarun);
        for cellid = datarun.cell_types{4}.cell_ids(1:20)
            figure();

            sanesubplot(2, 1, [1 1]);
            plot_rf(datarun, cellid);
            hold on;
            plot(cone_peaks(:,2), cone_peaks(:,1), 'yo');
            autozoom_to_fit(datarun, cellid, 5, 1, 1);

            sanesubplot(2, 1, [2 1]);
            imagesc(summed_rfs); colormap gray; axis equal tight ij
            hold on;
            plot(cone_peaks(:,2), cone_peaks(:,1), 'yo');
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
            plot(cone_peaks(:,1), cone_peaks(:,2), 'yo');
            autozoom_to_fit(datarun, cellid, 5, 1, 1);
            title(num2str(cellid));
        end


%% Manual checking and correction of cone locations
[new_locations, clicked_locations] = manual_cone_finding(datarun, {5}, cone_peaks, 'rad', 1);
% append the clicked (peak stixels) to cone_locations
cone_peaks = [cone_peaks; clicked_locations];   

% click to S cones from SBCs
gun = 3;
[new_locations, clicked_locations] = manual_cone_finding(datarun, {5}, new_locations, 'gun', gun, 'rad', 2);
        
% append the clicked (peak stixels) to cone_locations
cone_peaks = [cone_peaks; clicked_locations]; 

%% Optimize the cone center locations

% Find the center of mass for each cone, letting each cell that samples
% from the cone to contribute to this estimate
opt_cone_locations = optimize_cone_centers(datarun, cellspec, new_locations, all_sig_stix,...
                        'rad', gen_params.cone_location_optimization_roi);


%% Estimate the average cone RF size
cone_std = estimate_cone_size(datarun, cellspec, cone_peaks, opt_cone_locations, all_sig_stix,...
                        'rad', gen_params.cone_size_roi, 'upsampling', gen_params.cone_upsampling);


%% Make BW cone RFs

%intialize cone matrix
Wc = zeros([datarun.stimulus.field_width*datarun.stimulus.field_height], size(opt_cone_locations,1));

for cn = 1:size(opt_cone_locations,1)
    
    % set parameters to make cone rf
    cone_params.center_radius = cone_std;
    cone_params.center_scale = 1;
    cone_params.sparse = true;
    cone_params.x_size = datarun.stimulus.field_width;
    cone_params.y_size = datarun.stimulus.field_height;
    cone_params.effective_radius = ceil(3*cone_std);
    cone_params.center = [opt_cone_locations(cn,2), opt_cone_locations(cn,1)];
    
    temp_cone_rf = make_gaussian(cone_params);
    
    Wc(:,cn) = reshape(temp_cone_rf,[],1);
end

%% Estimate the cone RGB
cone_rgb = estimate_cone_rgb(datarun, cellspec, cone_peaks, opt_cone_locations, all_sig_stix, Wc);

%% If Color STAS then classify cones

if strcmpi(datarun.stimulus.independent, 't')
    classify_params.algorithm = 'k means';
    % classify cones
    [cone_types,likelihood_ratios,classify_extras] = ...
        classify_cones(cone_rgb, cone_rgb_expected(datarun),classify_params);
    nan_indices = isnan(cone_rgb(:,1))
    cone_types(nan_indices) = 'U';
else
    % chose a phony-baloney cone classification algorithm
    classify_params.algorithm = 'nearest line';
    cone_types = repmat('U', size(opt_cone_locations,1), 1);
    classify_extras.cone_rgb_observed = struct('U', [1 1 1]);
end


%% Identify cone sampling by each RGC

% free up memory prior to this computation
%clear all_sig_stixels cone_ideal_shapes

[cone_weights, regr_cell_ids,Wc]= extract_cone_weights(datarun, cellspec,...
    Wc, cone_types, classify_extras.cone_rgb_observed);


%% Fit center and surround gaussians to each RGC's cone RF
clear cone_info
cone_info.cone_weights = cone_weights;
cone_info.cone_centers = opt_cone_locations;
cone_info.cell_ids = regr_cell_ids;

rf_cone_fits = fit_cone_rfs(datarun, cellspec, 'cone_info', cone_info, 'fit_radius', 150, 'verbose', 1, 'foa_2d', []);


%% Save stuff out
addname = ['localmax-V2-st' num2str(stixel_threshold)];
savepath = sprintf('%s%s-%s/', single_cone_path, datarun.names.short_name, addname);
mkdir(savepath)

all_data.cone_weights = cone_weights;
all_data.cone_types = cone_types;
all_data.cone_centers = [opt_cone_locations(:,2), opt_cone_locations(:,1)];
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
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

save([savepath 'Wc'],'Wc')
save([savepath 'gen_params'],'gen_params')

%% Save some pictures 

% sensitivity surface and w/ cones superimposed
figure(1); clf
imagesc(summed_rfs); colormap gray;
axis off
print(1, [savepath,'sensitivity-surf.pdf'],'-dpdf')
figure(2); clf;
imagesc(summed_rfs); colormap gray; hold on
plot(opt_cone_locations(:,2), opt_cone_locations(:,1), 'g.');
axis off
hold off
print(2, [savepath,'sensitivity-surf-cones.pdf'],'-dpdf')

% plot the spiders for each of the major cell types
spider_threshold = 0.2; % percent of peak cone

%% Setup for Spider plots
%get size and color
y = datarun.stimulus.field_height;
x = datarun.stimulus.field_width;
rgb = [.5 .5 .5];
cone_halo_colors = [1 1 1]; 
cone_halo_colors = repmat(cone_halo_colors,4,1);

% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));

clear selection_params
selection_params.thresh = spider_threshold;
selection_params.radius = [0 inf];
selection_params.contiguity = true;
selection_params.polarity = 1;
halo_size = 6;
cone_size = 4;
line_max = 1.5;

%% on parasol
fig_num = 3;
spec = {1};
figure(fig_num);clf;image(plot_mat);axis image; hold on
[weights, center_selection, extras] = select_cone_weights(datarun, spec, selection_params);
summed_selection = sum(center_selection, 2);
sampled_indices = find(summed_selection > 0);
cone_center_roi = zeros(length(summed_selection),1);
cone_center_roi(sampled_indices) = 1;

plot_cell_sampling(datarun,spec, 'type', 'spider','fig_or_axes', fig_num, 'plot_cones', false,...
            'clear', false, 'label', false, 'thresh', spider_threshold,...
            'line_width',[realmin line_max], 'cell_colors', [1 1 1])
% plot cone halos
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',halo_size,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_halo_colors)
% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',cone_size,'clear',0,'bg_color',[])
title(['on-parasols: ',num2str(spider_threshold),'% of peak'])
print(fig_num, [savepath,'on-parasol-spiders.pdf'],'-dpdf')



%% off parsol 
fig_num = 4;
spec = {2};
figure(fig_num);clf;image(plot_mat);axis image; hold on
[weights, center_selection, extras] = select_cone_weights(datarun, spec, selection_params);
summed_selection = sum(center_selection, 2);
sampled_indices = find(summed_selection > 0);
cone_center_roi = zeros(length(summed_selection),1);
cone_center_roi(sampled_indices) = 1;

plot_cell_sampling(datarun,spec, 'type', 'spider','fig_or_axes', fig_num, 'plot_cones', false,...
            'clear', false, 'label', false, 'thresh', spider_threshold,...
            'line_width',[realmin line_max], 'cell_colors', [1 1 1])
% plot cone halos
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',halo_size,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_halo_colors)
% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',cone_size,'clear',0,'bg_color',[])
title(['off-parasols: ',num2str(spider_threshold),'% of peak'])
print(fig_num, [savepath,'off-parasol-spiders.pdf'],'-dpdf')


%% on midget
fig_num = 5;
spec = {3};
figure(fig_num);clf;image(plot_mat);axis image; hold on
[weights, center_selection, extras] = select_cone_weights(datarun, spec, selection_params);
summed_selection = sum(center_selection, 2);
sampled_indices = find(summed_selection > 0);
cone_center_roi = zeros(length(summed_selection),1);
cone_center_roi(sampled_indices) = 1;

plot_cell_sampling(datarun,spec, 'type', 'spider','fig_or_axes', fig_num, 'plot_cones', false,...
            'clear', false, 'label', false, 'thresh', spider_threshold,...
            'line_width',[realmin line_max], 'cell_colors', [1 1 1])
% plot cone halos
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',halo_size,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_halo_colors)
% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',cone_size,'clear',0,'bg_color',[])
title(['on-midgets: ',num2str(spider_threshold),'% of peak'])
print(fig_num, [savepath,'on-midget-spiders.pdf'],'-dpdf')


%% off midget
fig_num = 6;
spec = {4};
figure(fig_num);clf;image(plot_mat);axis image; hold on
[weights, center_selection, extras] = select_cone_weights(datarun, spec, selection_params);
summed_selection = sum(center_selection, 2);
sampled_indices = find(summed_selection > 0);
cone_center_roi = zeros(length(summed_selection),1);
cone_center_roi(sampled_indices) = 1;

plot_cell_sampling(datarun,spec, 'type', 'spider','fig_or_axes', fig_num, 'plot_cones', false,...
            'clear', false, 'label', false, 'thresh', spider_threshold,...
            'line_width',[realmin line_max], 'cell_colors', [1 1 1])
% plot cone halos
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',halo_size,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_halo_colors)
% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',cone_size,'clear',0,'bg_color',[])
title(['off-midgets: ',num2str(spider_threshold),'% of peak'])
print(fig_num, [savepath,'off-midget-spiders.pdf'],'-dpdf')


%% sbc 
sbc_thresh = 0.3;
fig_num = 7;
spec = {5};
sbc_selection_params.thresh = sbc_thresh;
sbc_selection_params.radius = [0 inf];
sbc_selection_params.contiguity = true;
sbc_selection_params.polarity = 1;
figure(fig_num);clf;image(plot_mat);axis image; hold on
[weights, center_selection, extras] = select_cone_weights(datarun, spec, selection_params);
summed_selection = sum(center_selection, 2);
sampled_indices = find(summed_selection > 0);
cone_center_roi = zeros(length(summed_selection),1);
cone_center_roi(sampled_indices) = 1;

plot_cell_sampling(datarun,spec, 'type', 'spider','fig_or_axes', fig_num, 'plot_cones', false,...
            'clear', false, 'label', false, 'thresh', sbc_thresh,...
            'line_width',[realmin line_max], 'cell_colors', [1 1 1], 'contiguity', false)
% plot cone halos
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',halo_size,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_halo_colors)
% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',cone_size,'clear',0,'bg_color',[])
title(['sbcs: ',num2str(spider_threshold),'% of peak'])
print(fig_num, [savepath,'sbc-spiders.pdf'],'-dpdf')



