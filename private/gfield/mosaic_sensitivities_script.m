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

plot_cell_sampling(datarun, {3})

datarun = select_cone_mosaic_roi(datarun, 'fig_or_axes', 12);

if 1
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


 % extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {1},...
                                            'thresh', 0.05,...
                                            'radius', [0 4], 'polarity', 0,...
                                            'contiguity', true,'scale', 3.0,...
                                            'remove_cones', 'U');   
        
connectivity = mosaic_weights .* selection;
num_cones = size(connectivity,1);

% normalize weights across rgcs by the peak cone
if 1
    norm_factors = max(connectivity, [], 1);
    connectivity = connectivity ./ repmat(norm_factors, num_cones, 1);
end



sensitivity_surface = sum(connectivity, 2) .* datarun.cones.roi;

sensitivity_image = Wc .* repmat(sensitivity_surface', size(Wc,1), 1);
sensitivity_image = sum(sensitivity_image,2);
figure
imagesc(norm_image(reshape(full(sensitivity_image), [320,320,3])))


figure(4)
sensitivity_surface(s_indices) = 0;
non_zero_inds = find(sensitivity_surface);
plot(sensitivity_surface(non_zero_inds))
coef_var = std(sensitivity_surface(non_zero_inds)) ./ mean(sensitivity_surface(non_zero_inds))


one_mosaic = sensitivity_surface(non_zero_inds);
mean(one_mosaic)
std(one_mosaic)


%%
 
% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {1,2,3,4},...
                                            'thresh', 0.05,...
                                            'radius', [0 4], 'polarity', 0,...
                                            'contiguity', false,'scale', 3.0,...
                                            'remove_cones', 'U');   
                
connectivity = mosaic_weights .* selection;

% normalize weights across rgcs by the peak cone
if 1
    norm_factors = max(connectivity, [], 1);
    connectivity = connectivity ./ repmat(norm_factors, num_cones, 1);
end


sensitivity_surface = sum(connectivity, 2) .* datarun.cones.roi;

sensitivity_image = Wc .* repmat(sensitivity_surface', size(Wc,1), 1);
sensitivity_image = sum(sensitivity_image,2);
figure
imagesc(norm_image(reshape(full(sensitivity_image), [320,320,3])))



figure
sensitivity_surface(s_indices) = 0;
non_zero_inds = find(sensitivity_surface);
plot(sensitivity_surface(non_zero_inds))
coef_var = std(sensitivity_surface(non_zero_inds)) ./ mean(sensitivity_surface(non_zero_inds))

remaining_mosaics = sensitivity_surface(non_zero_inds);
%%
figure
plot(one_mosaic, remaining_mosaics, 'k.')

figure
plot(one_mosaic)


figure
plot(remaining_mosaics)


one_mosaic_mean = mean(one_mosaic)
one_mosaic_std = std(one_mosaic)

temp_indices = find(one_mosaic < (one_mosaic_mean - (1 * one_mosaic_std)));

a = mean(one_mosaic);
b = mean(remaining_mosaics);

figure
plot(one_mosaic(temp_indices), remaining_mosaics(temp_indices), 'k.')
hold on
plot([a a], [50 300],'r')
plot([0 400], [b,b],'r')

remaining_comp_mean = mean(remaining_mosaics(temp_indices))
remaining_mean = mean(remaining_mosaics)


            
%% 
% sanity check on uniformity effect
temp_samples = normrnd(23.6, 38.1, 414,1);
temp_samples_two = normrnd(32.1, 33.7, 414,1);
temp_samples_three = normrnd(68.25, 71.88, 414,1);
temp_samples_four = normrnd(62.17, 60.63, 414,1);

% sanity check with normalization
temp_samples = normrnd(0.2022, 0.4416, 317,1);
temp_samples_two = normrnd(0.1808, 0.2487, 317,1);
temp_samples_three = normrnd(0.2244, 0.3479, 317,1);
temp_samples_four = normrnd(0.2452, 0.3390, 317,1);


summed_samples = temp_samples +temp_samples_two + temp_samples_three + temp_samples_four;

std(temp_samples) ./ mean(temp_samples)
std(temp_samples_two) ./ mean(temp_samples_two)
std(summed_samples) ./ mean(summed_samples)

% Indeed, expected CV falls as 1/sqrt(sample_number)
%%

datarun = fit_cone_rfs(datarun,{1,2,3,4});


cell_types = {4};
rgc_indices = get_cell_indices(datarun, cell_types);
num_rgcs = length(rgc_indices);
num_cones = length(datarun.cones.types);


% get weights from RF fits
new_mosaic_weights = zeros(num_cones, num_rgcs);


for rgc = 1:num_rgcs
    fit_params = [];
    fit_params.center = datarun.cones.rf_fits{rgc_indices(rgc)}.center;
    fit_params.center_scale = datarun.cones.rf_fits{rgc_indices(rgc)}.center_scale;
    fit_params.center_radius = datarun.cones.rf_fits{rgc_indices(rgc)}.center_radius;
    fit_params.surround_scale = datarun.cones.rf_fits{rgc_indices(rgc)}.surround_scale;
    fit_params.surround_radius = datarun.cones.rf_fits{rgc_indices(rgc)}.surround_radius;
    
    temp_weights = dog_fit(datarun.cones.centers, fit_params);

    new_mosaic_weights(:, rgc) = temp_weights;
end


connectivity = new_mosaic_weights;

% normalize weights across rgcs by the peak cone
if 0
    norm_factors = max(connectivity, [], 1);
    connectivity = connectivity ./ repmat(norm_factors, num_cones, 1);
end



sensitivity_surface = sum(connectivity, 2) .* datarun.cones.roi;

sensitivity_image = Wc .* repmat(sensitivity_surface', size(Wc,1), 1);
sensitivity_image = sum(sensitivity_image,2);
figure
imagesc(norm_image(reshape(full(sensitivity_image), [320,320,3])))


figure
sensitivity_surface(s_indices) = 0;
non_zero_inds = find(sensitivity_surface);
plot(sensitivity_surface(non_zero_inds))
coef_var = std(sensitivity_surface(non_zero_inds)) ./ mean(sensitivity_surface(non_zero_inds))





                                        
                                        


