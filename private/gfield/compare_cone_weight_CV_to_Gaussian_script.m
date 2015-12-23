% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';
cd /snle/lab/Experiments/Array/Shared/one/plantain/

% peach
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
path_and_name{1,2} = '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_20.00--standard';
cd /snle/lab/Experiments/Array/Shared/one/peach/

% apple
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data013/data013';
path_and_name{1,2} = '2010-03-05-2_rf-13-apple-gf-bayes-msf_25.00--standard';
cd /snle/lab/Experiments/Array/Shared/one/apple/

% kiwi
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{1,2} = '2008-05-13-3_data006_data006-bayes-msf_85.00--standard';
cd /snle/lab/Experiments/Array/Shared/one/kiwi/

% blueberry
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{1,2} = '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_20.00--standard';
cd /snle/lab/Experiments/Array/Shared/one/blueberry/


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

load Wc
%%

plot_cell_sampling(datarun, {4})

datarun = select_cone_mosaic_roi(datarun, 'fig_or_axes', 12);

if 0
    cd ~/Desktop/
    load cone_indices
    %cone_indices = cone_indices + 1;
end
temp_roi = zeros(length(datarun.cones.types),1);
temp_roi(cone_indices) =1;
datarun.cones.roi = temp_roi;

%%
% build sensitivity surface across mosaics

% find S cone indices to exclude from Wc
s_indices = find(datarun.cones.types == 'S');

cell_type = 3;
temp_cell_indices = get_cell_indices(datarun, {cell_type});

 % extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', 0.0,...
                                            'radius', [0 4], 'polarity', 0,...
                                            'contiguity', true,'scale', 3.0,...
                                            'remove_cones', 'U');   
        
num_cones = size(mosaic_weights,1);

% normalize weights across rgcs by the peak cone
if 1
    norm_factors = max(mosaic_weights, [], 1);
    mosaic_weights = mosaic_weights ./ repmat(norm_factors, num_cones, 1);
end

connectivity = mosaic_weights .* selection;

%collect the weights as a function from RF COM
figure(101); clf; hold on
for rgc = 1:length(temp_cell_indices)
    temp_dists = ipdm(datarun.stas.rf_coms{temp_cell_indices(rgc)}, datarun.cones.centers);
    temp_weights = mosaic_weights(:,rgc);
%    temp_weights = datarun.cones.weights(:,temp_cell_indices(rgc));
    plot(temp_dists, temp_weights, 'k.')
    temp_inds = find(temp_dists > 58 & temp_dists < 62);
    if rgc == 1
        distant_weights = temp_weights(temp_inds);
    else
        distant_weights = [distant_weights; temp_weights(temp_inds)];
    end
end
axis([0 80 -0.2 1.2])
sta_noise = std(distant_weights)

sensitivity_surface = sum(connectivity, 2) .* datarun.cones.roi;

% sensitivity_image = Wc .* repmat(sensitivity_surface', size(Wc,1), 1);
% sensitivity_image = sum(sensitivity_image,2);
% figure
% imagesc(norm_image(reshape(full(sensitivity_image), [320,320,3])))
% 
% 
figure(4)
sensitivity_surface(s_indices) = 0;
non_zero_inds = find(sensitivity_surface);
plot(sensitivity_surface(non_zero_inds))
coef_var = std(sensitivity_surface(non_zero_inds)) ./ mean(sensitivity_surface(non_zero_inds))


one_mosaic = sensitivity_surface(non_zero_inds);
mean(one_mosaic)
std(one_mosaic)

%% substitute Gaussian weights and recompute


temp_cell_indices = get_cell_indices(datarun, {cell_type});
Gauss_datarun = datarun;

for rgc = 1:length(temp_cell_indices)
    temp_fit = datarun.cones.rf_fits{temp_cell_indices(rgc)};
    rf_com = temp_fit.center;
    center_rad = temp_fit.center_radius;
    surround_rad = temp_fit.surround_radius;
    center_sc = temp_fit.center_scale;
    surround_sc = temp_fit.surround_scale;
    
    cone_dists = ipdm(rf_com, datarun.cones.centers);
    
    center_weights = center_sc * exp(-0.5.*(cone_dists.^2)./(center_rad.^2));
    surround_weights = surround_sc * exp(-0.5*(cone_dists.^2)./(surround_rad.^2));
    new_weights = center_weights - surround_weights;
    % add noise
    weight_noise = normrnd(repmat(0.0, size(new_weights)), repmat((max(new_weights) * sta_noise), size(new_weights)));
    Gauss_datarun.cones.weights(:,temp_cell_indices(rgc)) = (new_weights + weight_noise);
end
    

 % extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(Gauss_datarun, {cell_type},...
                                            'thresh', 0.0,...
                                            'radius', [0 4], 'polarity', 0,...
                                            'contiguity', true,'scale', 3.0,...
                                            'remove_cones', 'U');   

                                        
                                        
% normalize weights across rgcs by the peak cone
if 1
    norm_factors = max(mosaic_weights, [], 1);
    mosaic_weights = mosaic_weights ./ repmat(norm_factors, num_cones, 1);
end

connectivity = mosaic_weights .* selection;
num_cones = size(connectivity,1);

%collect the weights as a function from RF COM
figure(102); clf; hold on
for rgc = 1:length(temp_cell_indices)
    temp_dists = ipdm(datarun.stas.rf_coms{temp_cell_indices(rgc)}, Gauss_datarun.cones.centers);
    temp_weights = mosaic_weights(:,rgc);
    plot(temp_dists, temp_weights, 'k.')
    temp_inds = find(temp_dists > 58 & temp_dists < 62);
    if rgc == 1
        distant_weights = temp_weights(temp_inds);
    else
        distant_weights = [distant_weights; temp_weights(temp_inds)];
    end
end
axis([0 80 -0.2 1.2])
std(distant_weights)



sensitivity_surface = sum(connectivity, 2) .* Gauss_datarun.cones.roi;
% 
% sensitivity_image = Wc .* repmat(sensitivity_surface', size(Wc,1), 1);
% sensitivity_image = sum(sensitivity_image,2);
% figure
% imagesc(norm_image(reshape(full(sensitivity_image), [320,320,3])))


figure(4)
sensitivity_surface(s_indices) = 0;
non_zero_inds = find(sensitivity_surface);
plot(sensitivity_surface(non_zero_inds))
coef_var = std(sensitivity_surface(non_zero_inds)) ./ mean(sensitivity_surface(non_zero_inds))


one_mosaic = sensitivity_surface(non_zero_inds);
mean(one_mosaic)
std(one_mosaic)



    
    
    
    
    
    













