data_path = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
data_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';

cell_type = {1 2};

% load data
clear temp_datarun datarun
temp_datarun = load_data(data_path);
temp_datarun = load_params(temp_datarun,struct('verbose',1));  
temp_datarun = load_sta(temp_datarun, 'load_sta', cell_type);
datarun = temp_datarun;

datarun = get_sta_summaries(datarun, cell_type);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun);

temp_profiles = plot_rf_portraits(datarun, {2}, 'plot_radius', 40);



plot_cell_sampling(datarun, cell_type, 'label', true, 'cone_size', 14)

temp_profiles = plot_rf_portraits(datarun, [183 7696], 'plot_radius', 10,...
    'scale_factor', 10, 'cones', true, 'cone_size', 16, 'plot_color', []);

plot_cell_sampling(datarun, [183 7696], 'label', true)

print(2,'/snle/home/gfield/Desktop/profiles','-dpdf')

temp_index = get_cell_indicies(datarun, 183);
temp_rf = datarun.stas.rfs(datarun, 

figure(1)
x = [9 5 1];
explode = [0 0 1];
pie3(x, explode)
colormap([1 0 0; 0 1 0; 0 0 1])
print(1, '/snle/home/gfield/Desktop/center','-dpdf')

figure(2)
x = [29 33];
explode = [0 0];
pie3(x, explode)
colormap([1 0 0; 0 1 0])
print(2, '/snle/home/gfield/Desktop/surround','-dpdf')


temp_h = plot_cell_sampling(datarun, {3}, 'line_width', [0.1 2], 'cone_size', 5);
print(1, '/snle/home/gfield/Desktop/spider','-dpdf')


temp_params.criterion = 'snr';
temp_params.polarity = 1;
temp_params.thresh = 2.5;
temp_params.line_width = [0.1 2.0];
temp_params.cone_size = 5;

temp_h = plot_cell_sampling(datarun, {4});

%%%%%%%%%%%%%%%%%%%%%%%%
[286 6931 348 7607 901 1099 2821]

sigma_of_interest = 6;
temp_cell_index = get_cell_indices(datarun, 7607);
cell_fit = datarun.cones.rf_fits{temp_cell_index};
cell_center = cell_fit.center;
num_cones = length(datarun.cones.types);

cone_distances = zeros(1,num_cones);
for cone = 1:num_cones
    cone_center = datarun.cones.centers(cone,:);
    cone_distances(cone) = sqrt((cell_center(1) - cone_center(1)).^2 + (cell_center(2) - cone_center(2)).^2);
end

cones_of_interest = find(cone_distances < (cell_fit.center_radius .* sigma_of_interest));

[temp_max, temp_index] = max(datarun.cones.weights(:,temp_cell_index));

L_indices = find(datarun.cones.types(cones_of_interest) == 'L');
M_indices = find(datarun.cones.types(cones_of_interest) == 'M');
cell_weights = datarun.cones.weights(cones_of_interest,temp_cell_index);
L_weights = cell_weights(L_indices);
M_weights = cell_weights(M_indices);

net_L_weight = sum(L_weights)
net_M_weight = sum(M_weights)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma_of_interest = 4;
cell_indices = get_cell_indices(datarun, [286 6931 348 7607 901 1099 2821]);
num_cells = length(cell_indices);
num_cones = length(datarun.cones.types);
net_L_weights = zeros(1,num_cells);
net_M_weights = zeros(1,num_cells);
    
for cll = 1:num_cells
    cell_fit = datarun.cones.rf_fits{cell_indices(cll)};
    cell_center = cell_fit.center;
    cone_distances = zeros(1,num_cones);

    for cone = 1:num_cones
        cone_center = datarun.cones.centers(cone,:);
        cone_distances(cone) = sqrt((cell_center(1) - cone_center(1)).^2 + (cell_center(2) - cone_center(2)).^2);
    end

    cones_of_interest = find(cone_distances < (cell_fit.center_radius .* sigma_of_interest));
    length(cones_of_interest)
    
    L_indices = find(datarun.cones.types(cones_of_interest) == 'L');
    M_indices = find(datarun.cones.types(cones_of_interest) == 'M');
    cell_weights = datarun.cones.weights(cones_of_interest,temp_cell_index);
    L_weights = cell_weights(L_indices);
    M_weights = cell_weights(M_indices);

    temp_net_L_weight = sum(L_weights);
    temp_net_M_weight = sum(M_weights);
    
    net_L_weights(cll) = temp_net_L_weight;
    net_M_weights(cll) = temp_net_M_weight;

    %net_L_weights(cll) = temp_net_L_weight ./ (abs(temp_net_L_weight) + abs(temp_net_M_weight));
    %net_M_weights(cll) = temp_net_M_weight./ (abs(temp_net_L_weight) + abs(temp_net_M_weight));
end


plot(net_L_weights, net_M_weights, 'k.')


