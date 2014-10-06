%%
% MASTER
data_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data008/data008/data008';
obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/rf-8-gf/';

datarun{1} = load_data(data_path);
datarun{1} = load_params(datarun{1});
datarun{1} = load_neurons(datarun{1});
datarun{1} = load_ei(datarun{1}, 'all','array_type', 519);
datarun{1} = load_sta(datarun{1},'load_sta',[]);


datarun{1} = load_index(datarun{1});
datarun{1}.names.obvius_fit_path = obvius_fit_path;

datarun{1} = set_polarities(datarun{1});
datarun{1} = load_obvius_sta_fits(datarun{1});
datarun{1} = get_sta_fits_from_obvius(datarun{1}, {1,2,3,4});

% SLAVE
% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';
%path_and_name{1,2} = 'erroneous_normalization/2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard-old';

% load data
datarun{2} = load_data(path_and_name{1,1});
datarun{2} = load_params(datarun{2},struct('verbose',1));  
datarun{2} = load_sta(datarun{2},'load_sta',[]);
datarun{2} = load_index(datarun{2});
datarun{2} = load_ei(datarun{2}, 'all','array_type', 519);
datarun{2} = load_sta(datarun{2},'load_sta',[]);

datarun{2} = get_sta_summaries(datarun{2}, {1,2,3,4}, 'keep_stas', false);

%datarun{2} = set_polarities(datarun{2});
datarun{2} = import_single_cone_data(datarun{2}, path_and_name{1,2});    
datarun{2}.cones.mosaic = make_mosaic_struct(datarun{2}.cones.centers);

mapped_cell_list = map_ei(datarun{1}, datarun{2}, 'master_cell_type', {3}, 'slave_cell_type', {3});

counter = 0;
for cc = 1:length(mapped_cell_list)
    if ~isempty(mapped_cell_list{cc})
        counter = counter + 1;
        master_cell_indices(counter) = cc;
        slave_cell_ids(counter) = mapped_cell_list{cc};
    end
end

slave_cell_indices = get_cell_indices(datarun{2}, slave_cell_ids);

% save_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data008/data008/';
% map_ei_classification_txt(datarun{1}, datarun{2},'master_cell_type',{3}, 'slave_cell_type', {3},...
%                         'classification_path', save_path, 'corr_threshold', 0.95);
% 
%                     
% 
% datarun{2}.names.map_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data008/data008/ei_map-from-.txt';
% datarun = load_map(datarun);

%%






cell_types = {3};
pixel_rescale = 5; % scale factor between 1x1 and coarse spatial STA run

num_rgcs = length(master_cell_indices);

cone_locations = datarun{2}.cones.centers;
new_datarun = datarun{2};
num_cones = size(cone_locations,1);
new_weights = zeros(num_cones,num_rgcs);

for rgc = 1:num_rgcs
    % extract and organized fit information
    the_fit = datarun{1}.stas.fits{master_cell_indices(rgc)};
    fit_params = [];
    fit_params.center = datarun{2}.stas.rf_coms{slave_cell_indices(rgc)};
    fit_params.center_scale = 1;
    fit_params.center_sd = the_fit.sd .* pixel_rescale;
    fit_params.angle = the_fit.angle;
    fit_params.surround_sd_scale = the_fit.surround_sd_scale;
    fit_params.surround_scale = the_fit.surround_scale;
    
    temp_weights = two_d_dog_fit(cone_locations, fit_params);
    new_weights(:,rgc) = temp_weights;
    new_datarun.cones.weights(:,slave_cell_indices(rgc)) = temp_weights;
    %new_datarun.stas.rf_coms{slave_cell_indices(rgc)} = fit_params.center;
    
%    recon_rf = Wc .* repmat(temp_weights', size(Wc,1), 1); 
%    recon_rf = sum(recon_rf,2);
%    figure
%    imagesc(norm_image(reshape(full(recon_rf), [320,320,3])))

    
end

    
figure_num = 11;

cell_spec = {3};
temp_polarity = 1;

% plot background

% get size and color
y = new_datarun.stimulus.field_height;
x = new_datarun.stimulus.field_width;
rgb = [.35 .35 .35];

% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));

% plot it
figure(figure_num);clf;image(plot_mat);axis image; hold on


clear selection_params
selection_params.thresh = 0.05;
selection_params.radius = [0 8];
selection_params.contiguity = true;
selection_params.polarity = 1;
[weights, center_selection, extras] = select_cone_weights(new_datarun, slave_cell_ids, selection_params);

summed_selection = sum(center_selection, 2);
sampled_indices = find(summed_selection > 0);
cone_center_roi = zeros(length(summed_selection),1);
cone_center_roi(sampled_indices) = 1;

% plot spiders
plot_cell_sampling(new_datarun,slave_cell_ids,'type','spider','fig_or_axes',figure_num,'clear',0,'plot_cones',0,...
   'line_width',[realmin 2.0],'cell_colors',...
    [0.5 0.5 0.5] + (0.5*temp_polarity),'plot_radius',[0 6],...
    'thresh', selection_params.thresh,'contiguity', selection_params.contiguity, 'polarity', 0)


% plot center cone halos
cone_colors = [0.5 0.5 0.5] + (0.5*temp_polarity); 
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(new_datarun,'fig_or_axes',figure_num,'cone_size',4,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(new_datarun,'fig_or_axes',figure_num,'cone_size',8,'clear',0,'bg_color',[])


%%

cone_dists = zeros(num_cones, num_rgcs);
for rgc = 1:num_rgcs
    temp_com = datarun{2}.stas.rf_coms{slave_cell_indices(rgc)};
    temp_dists = ipdm(cone_locations, temp_com);
    cone_dists(:,rgc) = temp_dists;
end

cone_dists = reshape(cone_dists, [],1);
profile_weights = reshape(new_weights, [],1);
figure
plot(cone_dists, profile_weights, 'k.');
axis([0 50 -0.1 1])


 % extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun{2}, slave_cell_ids,...
                                            'thresh', 0.00,...
                                            'radius', [0 inf], 'polarity', 0,...
                                            'contiguity', false,'scale', 3.0,...
                                            'remove_cones', 'U');   


connectivity = mosaic_weights .* selection;                                        
% normalize cone weights
for RGC = 1:num_rgcs
    connectivity(:,RGC) = connectivity(:,RGC) ./ datarun{2}.cones.rf_fits{slave_cell_indices(RGC)}.center_scale;
end                                         
            
observed_weights = reshape(connectivity, [],1);
%non_zero_indices = find(observed_weights ~= 0);

figure
plot(cone_dists, observed_weights, 'k.');
axis([0 20 -0.1 1])


%%
% fit the single cone RFs
%temp_datarun = fit_cone_rfs(datarun{2},{3});


for rgc = 1:num_rgcs
    temp_com = datarun{2}.stas.rf_coms{slave_cell_indices(rgc)};
    temp_dists = ipdm(cone_locations, temp_com);
    temp_weights = datarun{2}.cones.weights(:,slave_cell_indices(rgc));
    
    plot(temp_dists, temp_weights, 'k.')
    pause(0.5)
    hold on

end

figure
for rgc = 1:num_rgcs
    temp_weights = dog_fit(cone_locations, datarun{2}.cones.rf_fits{slave_cell_indices(rgc)});
    temp_dists = ipdm(cone_locations, datarun{2}.cones.rf_fits{slave_cell_indices(rgc)}.center);
    
    plot(temp_dists, temp_weights, 'k.')
    pause(0.1)
    hold on
end






