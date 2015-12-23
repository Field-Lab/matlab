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

% make blurred cone mosaic image
blurred_cone_mosaic = compute_blurred_cone_mosaic(datarun, 'verbose', true);

% plot come mosaic
background = zeros(size(blurred_cone_mosaic)) + 0.2;
figure(1); clf;
sat_cone_mosaic = (blurred_cone_mosaic + background) * 1.1;
sat_cone_mosaic(sat_cone_mosaic > 1) = 1.0;
image(sat_cone_mosaic)
axis image; axis off; hold on
title('original cone mosaic')

%%

% the 1-sigma fits over the top of this
cell_type = 4;
cell_indices = get_cell_indices(datarun, {cell_type});
num_rgcs = length(cell_indices);

% get connectivity and purities
% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0,...
                                            'remove_cones', 'U');   
connectivity = mosaic_weights .* selection; 
purity_indices = compute_opponency_index(connectivity, datarun.cones.types);
mean_purity = mean(purity_indices);
relative_purities = purity_indices;

% set up colormap for rf fits
rg_map = zeros(num_rgcs,3);
ln_widths = zeros(num_rgcs,1);
for rgc = 1:num_rgcs
    if relative_purities(rgc) < -0.8
        rg_map(rgc,:) = [1 1 1];
        ln_widths(rgc) = 1.5;
    elseif relative_purities(rgc) < -0.6 && relative_purities(rgc) >= -0.8
        rg_map(rgc,:) = [0.9 0.9 0.9];
        ln_widths(rgc) = 1.5;
    elseif relative_purities(rgc) < -0.4 && relative_purities(rgc) >= -0.6
        rg_map(rgc,:) = [0.8 0.8 0.8];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) < -0.2 && relative_purities(rgc) >= -0.4
        rg_map(rgc,:) = [0.7 0.7 0.7];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) <= 0.2 && relative_purities(rgc) >= -0.2
        rg_map(rgc,:) = [0.5 0.5 0.5];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) > 0.2 && relative_purities(rgc) <= 0.4
        rg_map(rgc,:) = [0.5 0.5 0.5];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) > 0.4 && relative_purities(rgc) <= 0.6
        rg_map(rgc,:) = [0.2 0.2 0.2];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) > 0.6 && relative_purities(rgc) <= 0.8
        rg_map(rgc,:) = [0.1 0.1 0.1];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) > 0.8
        rg_map(rgc,:) = [0 0 0];
        ln_widths(rgc) = 1.5;
    end
end

num_rgcs = length(cell_indices);

for rgc = 1:num_rgcs
    
        % get the fit
        the_fit = datarun.cones.rf_fits{cell_indices(rgc)};

        % skip if doesn't exist
        if isempty(the_fit);continue;end

        % get center
        ctr = the_fit.center;

        % get center radius
        rad = the_fit.center_radius;
        rad = [rad rad];
        
        fit_angle = 0;
        
        % get points of an ellipse with these parameters
        [X,Y] = drawEllipse([ctr rad fit_angle]);
        
        % if any are NaN, skip
        if any(isnan([X Y]));continue;end
        
        % transform to desired coordinates
 %       [X, Y] = tformfwd(coord_tform, X, Y);

        % plot the points
        plot(X,Y,'Color',rg_map(rgc,:), 'LineWidth', ln_widths(rgc))
end

%% Shuffle the cone IDs
L_indices = find(datarun.cones.types == 'L');
M_indices = find(datarun.cones.types == 'M');
S_indices = find(datarun.cones.types == 'S');

LM_indices = [L_indices; M_indices];
num_LM = length(LM_indices);
shuffled_LM = randperm(num_LM);

new_LM_indices = LM_indices(shuffled_LM);
cone_types = datarun.cones.types;

num_cones = length(datarun.cones.types);
new_cone_types(1:num_cones) = 'S';

new_cone_types(new_LM_indices) = cone_types(LM_indices);
new_cone_types = char(new_cone_types);

new_datarun = datarun;
new_datarun.cones.types = new_cone_types;

%% load simulated cone mosaics.
clumped_flag = true;
cone_mosaic = 1;

%clear cleaned_cone_types temp_types cone_types
suffix_string = num2str(11110+cone_mosaic);
if clumped_flag
    str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations/clumped/',path_and_name{1,2},'-',suffix_string,'.txt'];
else
    str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations/permuted/',path_and_name{1,2},'-',suffix_string,'.txt'];  
end
cones_file = dlmread([str2load],'\t');
cone_types = cones_file(:,4);
%S_indices = find(cone_types == 3); % s cones
%LM_cone_indices = setdiff([1:length(cone_types)], S_indices);    

%temp_types = cone_types(keep_cone_indices);
clumped_cone_types(cone_types == 1,1) = 'L';
clumped_cone_types(cone_types == 2,1) = 'M';
clumped_cone_types(cone_types == 3,1) = 'S';

permuted_cone_PIs = compute_opponency_index(connectivity, clumped_cone_types);
permuted_cone_PI_SD(cone_mosaic) = std(permuted_cone_PIs)

new_datarun = datarun;
new_datarun.cones.types = clumped_cone_types;

%%
% make blurred cone mosaic image
clumped_blurred_cone_mosaic = compute_blurred_cone_mosaic(new_datarun, 'verbose', true);

% plot come mosaic
background = zeros(size(clumped_blurred_cone_mosaic)) + 0.2;
figure(2); clf;
sat_cone_mosaic = (clumped_blurred_cone_mosaic + background) * 1.1;
sat_cone_mosaic(sat_cone_mosaic > 1) = 1.0;
image(sat_cone_mosaic)
axis image; axis off; hold on
title('clumped cone mosaic')

%%

% get connectivity and purities
% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(new_datarun, {cell_type},...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0,...
                                            'remove_cones', 'U');   
connectivity = mosaic_weights .* selection; 
purity_indices = compute_opponency_index(connectivity, new_datarun.cones.types);
mean_purity = mean(purity_indices);
relative_purities = purity_indices;

% set up colormap for rf fits
rg_map = zeros(num_rgcs,3);
ln_widths = zeros(num_rgcs,1);
for rgc = 1:num_rgcs
    if relative_purities(rgc) < -0.8
        rg_map(rgc,:) = [1 1 1];
        ln_widths(rgc) = 1.5;
    elseif relative_purities(rgc) < -0.6 && relative_purities(rgc) >= -0.8
        rg_map(rgc,:) = [0.9 0.9 0.9];
        ln_widths(rgc) = 1.5;
    elseif relative_purities(rgc) < -0.4 && relative_purities(rgc) >= -0.6
        rg_map(rgc,:) = [0.8 0.8 0.8];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) < -0.2 && relative_purities(rgc) >= -0.4
        rg_map(rgc,:) = [0.7 0.7 0.7];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) <= 0.2 && relative_purities(rgc) >= -0.2
        rg_map(rgc,:) = [0.5 0.5 0.5];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) > 0.2 && relative_purities(rgc) <= 0.4
        rg_map(rgc,:) = [0.5 0.5 0.5];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) > 0.4 && relative_purities(rgc) <= 0.6
        rg_map(rgc,:) = [0.2 0.2 0.2];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) > 0.6 && relative_purities(rgc) <= 0.8
        rg_map(rgc,:) = [0.1 0.1 0.1];
        ln_widths(rgc) = 1;
    elseif relative_purities(rgc) > 0.8
        rg_map(rgc,:) = [0 0 0];
        ln_widths(rgc) = 1.5;
    end
end


num_rgcs = length(cell_indices);

for rgc = 1:num_rgcs
    
        % get the fit
        the_fit = datarun.cones.rf_fits{cell_indices(rgc)};

        % skip if doesn't exist
        if isempty(the_fit);continue;end

        % get center
        ctr = the_fit.center;

        % get center radius
        rad = the_fit.center_radius;
        rad = [rad rad];
        
        fit_angle = 0;
        
        % get points of an ellipse with these parameters
        [X,Y] = drawEllipse([ctr rad fit_angle]);
        
        % if any are NaN, skip
        if any(isnan([X Y]));continue;end
        
        % transform to desired coordinates
 %       [X, Y] = tformfwd(coord_tform, X, Y);

        % plot the points
        plot(X,Y,'Color',rg_map(rgc,:), 'LineWidth', ln_widths(rgc))
end



%% load simulated cone mosaics.
clumped_flag = false;
cone_mosaic = 1;

%clear cleaned_cone_types temp_types cone_types
suffix_string = num2str(11110+cone_mosaic);
if clumped_flag
    str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations/clumped/',path_and_name{1,2},'-',suffix_string,'.txt'];
else
    str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations/permuted/',path_and_name{1,2},'-',suffix_string,'.txt'];  
end
cones_file = dlmread([str2load],'\t');
cone_types = cones_file(:,4);
%S_indices = find(cone_types == 3); % s cones
%LM_cone_indices = setdiff([1:length(cone_types)], S_indices);    

%temp_types = cone_types(keep_cone_indices);
clumped_cone_types(cone_types == 1,1) = 'L';
clumped_cone_types(cone_types == 2,1) = 'M';
clumped_cone_types(cone_types == 3,1) = 'S';

permuted_cone_PIs = compute_opponency_index(connectivity, clumped_cone_types);
permuted_cone_PI_SD(cone_mosaic) = std(permuted_cone_PIs)

new_datarun = datarun;
new_datarun.cones.types = clumped_cone_types;

%%
% make blurred cone mosaic image
random_blurred_cone_mosaic = compute_blurred_cone_mosaic(new_datarun, 'verbose', true);

% plot come mosaic
background = zeros(size(random_blurred_cone_mosaic)) + 0.2;
figure(3); clf;
sat_cone_mosaic = (random_blurred_cone_mosaic + background) * 1.1;
sat_cone_mosaic(sat_cone_mosaic > 1) = 1.0;
image(sat_cone_mosaic)
axis image; axis off; hold on
title('clumped cone mosaic')







