% peach
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
path_and_name{1,2} = 'peach';

% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain';
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data008/data008/data008';
path_and_name{1,2} = 'plantain';

% blueberry
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{1,2} = 'blueberry';



% load data
datarun = load_data(path_and_name{1,1});
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);

obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/rf-8-gf/stas-008/';
datarun.names.obvius_fit_path = obvius_fit_path;

datarun = load_obvius_sta_fits(datarun);
datarun = get_sta_fits_from_obvius(datarun, {1,2,3,4,5});


% get sta summaries 
datarun = get_sta_summaries(datarun, {3, 4}, 'keep_stas', false);

% portraits
cell_type = 5;
plot_rf_portraits(datarun, {cell_type}, 'plot_radius', 38);


% import cone information
%datarun = import_single_cone_data(datarun, path_and_name{1, 2});




%%%%%%%%%%%%%%%%%%
% Peach: cell list
on_parasol_one = 4025; on_parasol_two = 3856;
off_parasol_one = 4277; off_parasol_two = 3334;
on_midget_one = 3811; on_midget_two = 1516;
off_midget_one = 4022; off_midget_two = 3603;

%%%%%%%%%%%%%%%%%%
% Plantain: cell list
% data003 (1x1)
on_parasol_one = 1445; on_parasol_two = 1771;
off_parasol_one = 887; off_parasol_two = 526;
on_midget_one = 2117; on_midget_two = 3033;  %1741
off_midget_one = 976; off_midget_two = 1426;
sbc_one = 31; sbc_two = 603;

% data008 (5x5)
on_parasol_one = 1442; on_parasol_two = 1503;
off_parasol_one = 661; off_parasol_two = 1038;
on_midget_one = 2116; on_midget_two = 3031;  %1741
off_midget_one = 977; off_midget_two = 1487;
sbc_one = 7711; sbc_two = 602;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


window_size = 30; % 1x1
window_size = 6; % 5x5
% set default fit location
datarun.default_sta_fits = 'vision'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% on parasol rfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% on parasol mosaic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot_rf_fit(datarun, {1}, 'fits_to_use', 'obvius','fig_or_axes', 1)
axis([5 45 10 50])
axis square
print(1, '~/Desktop/on-parasol-mosaic','-deps')

% example cell 1
cell_one = on_parasol_one;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 1, 'scale', 8)

%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])

print(1, '~/Desktop/on-parasol-one','-dpdf')

% example cell 2
cell_two = on_parasol_two;
cell_two_index = get_cell_indices(datarun, cell_two);
plot_rf(datarun, cell_two, 'foa', 2, 'scale', 8)

%get COM
temp_COM = datarun.stas.rf_coms{cell_two_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])


print(2, '~/Desktop/on-parasol-two','-dpdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% off parasol rfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% off parasol mosaic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot_rf_fit(datarun, {2}, 'fits_to_use', 'obvius','fig_or_axes', 1)
axis([5 45 10 50])
axis square
print(1, '~/Desktop/off-parasol-mosaic','-deps')

% example cell 1
cell_one = off_parasol_one;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 2, 'scale', 8, 'polarity', true)

%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])

print(2, '~/Desktop/off-parasol-one','-dpdf')

% example cell 2
cell_two = off_parasol_two;
cell_two_index = get_cell_indices(datarun, cell_two);
plot_rf(datarun, cell_two, 'foa', 3, 'scale', 8, 'polarity', true)

%get COM
temp_COM = datarun.stas.rf_coms{cell_two_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])


print(3, '~/Desktop/off-parasol-two','-dpdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% on midget rfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% on midget mosaic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot_rf_fit(datarun, {3}, 'fits_to_use', 'obvius','fig_or_axes', 1,...
                'edge_color', 'w', 'fill', true, 'fill_color', [0.75 0.75 0.75])
axis([5 45 10 50])
axis square
print(1, '~/Desktop/on-midget-mosaic','-deps')

% example cell 1
cell_one = on_midget_one;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 2, 'scale', 8, 'polarity', true)

%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])

print(2, '~/Desktop/on_midget-one','-dpdf')

% example cell 2
cell_two = on_midget_two;
cell_two_index = get_cell_indices(datarun, cell_two);
plot_rf(datarun, cell_two, 'foa', 3, 'scale', 8, 'polarity', true)

%get COM
temp_COM = datarun.stas.rf_coms{cell_two_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])


print(3, '~/Desktop/on-midget-two','-dpdf')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% off midget rfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% off midget mosaic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot_rf_fit(datarun, {4}, 'fits_to_use', 'obvius','fig_or_axes', 1)
axis([5 45 10 50])
axis square
print(1, '~/Desktop/off-midget-mosaic','-dpdf')

% example cell 1
cell_one = off_midget_one;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 2, 'scale', 8, 'polarity', true)

%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])

print(2, '~/Desktop/off_midget-one','-dpdf')

% example cell 2
cell_two = off_midget_two;
cell_two_index = get_cell_indices(datarun, cell_two);
plot_rf(datarun, cell_two, 'foa', 3, 'scale', 8, 'polarity', true)

%get COM
temp_COM = datarun.stas.rf_coms{cell_two_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])


print(3, '~/Desktop/off-midget-two','-dpdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBC rfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% portraits
cell_type = 5;
plot_rf_portraits(datarun, {cell_type}, 'plot_radius', 45, 'plot_color', 3);

window_size = 40;

figure(1)
plot_rf_fit(datarun, {5})
axis([5 45 10 50])
axis square
print(1, '~/Desktop/sbc-mosaic','-dpdf')

% black and white of S cones
plot_rf_portraits(datarun, sbc_one, 'plot_radius', (8), 'plot_color', 3, 'scale_factor', 8)
plot_rf_portraits(datarun, sbc_two, 'plot_radius', (8), 'plot_color', 3, 'scale_factor', 8)

color_transform = [1 0 0; 0 1 0; 0.25 0.5 1];

% example cell 1
cell_one = sbc_one;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 2, 'scale', 8, 'com', false, 'color_transform', color_transform)

%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])

print(2, '~/Desktop/sbc-one','-dpdf')

% example cell 2
cell_two = sbc_two;
cell_two_index = get_cell_indices(datarun, cell_two);
plot_rf(datarun, cell_two, 'foa', 3, 'scale', 8, 'color_transform', color_transform)

%get COM
temp_COM = datarun.stas.rf_coms{cell_two_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])


print(3, '~/Desktop/sbc-two','-dpdf')


