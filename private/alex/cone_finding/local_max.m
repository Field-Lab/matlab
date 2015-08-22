%% LOCAL MAX CONE FINDING
% This code is an attempt to simpify and streamline earlier versions of
% local-max cone finding.  
% :GDF 2013-02

%% BW data only

% for SBC, fits do not work. Increase scale or do manual cone adding
path2load = '/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data008-from-d08_11/data008-from-d08_11';

frame = 4;
cell_specification = {5,3,4};
scale = 3;
field_size = [300 300];

%% load data
% main datarun
datarun = load_data(path2load);
datarun = load_params(datarun);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

[cell_ids, cell_type] = get_cell_indices(datarun, cell_specification);

ncells = length(datarun.cell_ids);

%% construct mosaics from fits

coord_tform = coordinate_transform(datarun,'sta');

ctr=nan(ncells,2);
rad=ctr;
fit_angle=nan(ncells,1);
for cell_index = 1:ncells
    the_fit = datarun.stas.fits{cell_index};
    if ~isempty(the_fit)
        ctr(cell_index,:) = the_fit.mean;
        rad(cell_index,:) = the_fit.sd;
        fit_angle(cell_index)=the_fit.angle;
    end
end

all_sta = zeros([field_size, length(cell2mat(cell_ids))]);
all_marks = all_sta;
all_rfs = zeros([field_size,length(cell_specification)]);
figure
cnt = 1;
for i=1:length(cell_specification)
    
    summed_rfs = 0;
    pol = datarun.stas.polarities{cell_ids{i}(1)};
    subplot(2,3,i)
    colormap gray
    
    for rgc=cell_ids{i}
        sta = datarun.stas.stas{rgc}(:,:,1,frame);
        all_sta(:,:,cnt) = sta;
        
        if ~isnan(rad(rgc,1))
            [X, Y] = drawEllipse([ctr(rgc,:) rad(rgc,:)*scale fit_angle(rgc)]);
            [X, Y] = tformfwd(coord_tform, X, Y);
            [inPoints] = polygrid( round(X), round(Y), 1);
            new_fit.ctr(cnt,:) = ctr(rgc,:);
            new_fit.rad(cnt,:) = rad(rgc,:)*scale;
            new_fit.fit_angle(cnt) = fit_angle(rgc);
            
            % cut borders if fit is beyond stim field
            inPoints(inPoints<1)=1;
            inPoints(inPoints(:,1)>field_size(1),1)=field_size(1);
            inPoints(inPoints(:,2)>field_size(2),2)=field_size(2);
            
            tmp = zeros(size(sta));
            for j=1:length(inPoints)
                tmp(inPoints(j,2),inPoints(j,1)) = pol*sta(inPoints(j,2),inPoints(j,1));
            end
            tmp = tmp/max(tmp(:));
            tmp(tmp<0) = 0; % cut negative pixels
            summed_rfs = summed_rfs + tmp;
            
            all_marks(:,:,cnt) = tmp;            
            
            imagesc(summed_rfs)
            title(int2str(i))
            drawnow
        end
        cnt = cnt+1;        
    end
    
    all_rfs(:,:,i) = summed_rfs/max(summed_rfs(:));
end
all_marks(all_marks>0) = 1;

all_sta = imresize(all_sta,2);
all_marks = imresize(all_marks,2);


figure
imagesc(all_rfs)

tmp = all_rfs(:,:,2)+all_rfs(:,:,3);
figure
colormap gray
imagesc(tmp)

comb = zeros(300,300,3);
comb(:,:,2) = tmp/max(tmp(:));
comb(:,:,3) = all_rfs(:,:,1);
figure
imagesc(imresize(comb,2, 'bicubic'))


final_sum = sum(imresize(comb(:,:,1:2),2, 'bicubic'),3);
figure
colormap gray
imagesc(final_sum)
hold on
cone_peaks = find_local_maxima(final_sum, 'radius', 4, 'thresh', 0.1, 'return', 'indices');
plot(cone_peaks(:,2), cone_peaks(:,1), 'xb')
set(gca, 'dataaspectratio',[1 1 1])


% make stimulation regions
diam = 3;
for i=1:length(cone_peaks)
    [X, Y] = drawEllipse([cone_peaks(i,:) [diam diam] pi/4]);
    [X, Y] = tformfwd(coord_tform, X, Y);
%     [inPoints] = polygrid( round(X), round(Y), 1);
%     plot(inPoints(:,2), inPoints(:,1), 'x')
    plot(round(Y),round(X),'r')
end

voronoi(cone_peaks(:,2), cone_peaks(:,1));
for i=1:40:600-40
    line([1 600], [i, i], 'color', 'g')
    line( [i, i],[1 600], 'color', 'g')
end


%% inspect
base_diam = 3;
cone_peaks = find_local_maxima(final_sum, 'radius', 4, 'thresh', 0.1, 'return', 'indices');


close all
cone_peaks1= inspect_cone_finding(all_sta, all_marks, final_sum, cone_peaks, base_diam, new_fit, coord_tform);


figure

for i=1:40:600-40
    for j=1:40:600-40
        
        % 
        
    end    
end






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



