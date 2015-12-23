% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain';

% apricot
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005-s3600-s7200/data005/data005';
path_and_name{1,2} = 'apricot';

% peach
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
path_and_name{1,2} = 'peach';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);


blur_param = 1;


%%

blurred_cone_mosaic = compute_blurred_cone_mosaic(datarun, 'verbose', true, 'blur_factor', blur_param);

background = zeros(size(blurred_cone_mosaic)) + 0.2;
figure; clf;
sat_cone_mosaic = (blurred_cone_mosaic + background) * 1.1;
sat_cone_mosaic(sat_cone_mosaic > 1) = 1.0;
image(sat_cone_mosaic)
axis image; axis off; hold on
title('original cone mosaic')


%% extra clumped
suffix_string = num2str(11111);
str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations/extra-clumped/',path_and_name{1,2},'-',suffix_string,'.txt'];
cones_file = dlmread([str2load],'\t');
cone_types = cones_file(:,4);

cleaned_cone_types(cone_types == 1,1) = 'L';
cleaned_cone_types(cone_types == 2,1) = 'M';
cleaned_cone_types(cone_types == 3,1) = 'S';
new_datarun = datarun;

new_datarun.cones.types = cleaned_cone_types;


%%

blurred_cone_mosaic_sim = compute_blurred_cone_mosaic(new_datarun, 'verbose', true, 'blur_factor', blur_param);

background = zeros(size(blurred_cone_mosaic_sim)) + 0.2;
figure; clf;
sat_cone_mosaic = (blurred_cone_mosaic_sim + background) * 1.1;
sat_cone_mosaic(sat_cone_mosaic > 1) = 1.0;
image(sat_cone_mosaic)
axis image; axis off; hold on
title('extra-clumped')


%% random
suffix_string = num2str(11111);
str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations/permuted/',path_and_name{1,2},'-',suffix_string,'.txt'];
cones_file = dlmread([str2load],'\t');
cone_types = cones_file(:,4);

cleaned_cone_types(cone_types == 1,1) = 'L';
cleaned_cone_types(cone_types == 2,1) = 'M';
cleaned_cone_types(cone_types == 3,1) = 'S';
random_datarun = datarun;

random_datarun.cones.types = cleaned_cone_types;


%%


blurred_cone_mosaic_sim = compute_blurred_cone_mosaic(random_datarun, 'verbose', true, 'blur_factor', blur_param);

background = zeros(size(blurred_cone_mosaic_sim)) + 0.2;
figure; clf;
sat_cone_mosaic = (blurred_cone_mosaic_sim + background) * 1.1;
sat_cone_mosaic(sat_cone_mosaic > 1) = 1.0;
image(sat_cone_mosaic)
axis image; axis off; hold on
title('random-permuted')

