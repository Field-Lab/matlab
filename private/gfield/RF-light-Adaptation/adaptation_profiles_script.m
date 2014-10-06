%%
data_paths{1} = '2007-03-08-4/data003/data003'; % ndf 2.0
data_paths{2} = '2007-03-08-4/data005/data005'; % ndf 0.0

%%
data_paths{1} = '2007-03-08-3/data011-gdf/data011'; % ndf 3.3  % good for on and off parasols
data_paths{2} = '2007-03-08-3/data012-gdf/data012'; % ndf 0.0

%%
data_paths{1} = '2007-03-08-1/data002-gdf/data002'; % ndf 2.6  % good for on and off parasols
data_paths{2} = '2007-03-08-1/data004-gdf/data004'; % ndf 0.0

%%
data_paths{1} = '2007-03-01-2/data002-gdf/data002'; % ndf 3.9  % good for on and off parasols
data_paths{2} = '2007-03-01-2/data005-gdf/data005'; % ndf 0.0

%%
data_paths{1} = '2007-03-02-0/data001/data001'; % ndf 3.9  % good for on parasols and midgets  % Maybe usable
data_paths{2} = '2007-03-02-0/data015/data015'; % ndf 0.0

%%
data_paths{1} = '2007-01-23-5/data001/data001'; % ndf 4.0  % good for on parasols and midgets  % Maybe usable
data_paths{2} = '2007-01-23-5/data010/data010'; % ndf 0.0

%%
data_paths{1} = '2007-03-27-2/data001-gdf/data001'; % ndf 3.9  % good for on parasols and midgets  % Maybe usable
data_paths{2} = '2007-03-27-2/data009/data009-nwpca/data009'; % ndf 0.0

%%
data_paths{1} = '2009-04-13-2/data006/data006'; % ndf 4.0  % OK for on parasols % Maybe usable
data_paths{2} = '2009-04-13-2/data011/data011'; % ndf 0.0

%%
data_paths{1} = '2008-04-22-3/data002/data002'; % ndf 4.3  % parasol SNR is low
%data_paths{2} = '2008-04-22-3/data006/data006'; % ndf 0.0

%%
data_paths{1} = '2008-04-08-2/data000-mg/data000'; % ndf 4.0  % parasols are good, midgets may be too coarse
data_paths{1} = '2008-04-08-2/data002/data002'; % ndf 4.0  % parasols are too incomplete, maybe use data003 for parasols
%data_paths{2} = '2008-04-08-2/data004/data004'; % ndf 0.0

%%
data_paths{1} = '2007-03-27-1/data007/data007'; % ndf 3.6  % Very Good
data_paths{2} = '2007-03-27-1/data010/data010'; % ndf 0.6

base_fit_path = '/snle/lab/Experiments/Array/Analysis/2007-03-27-1/';
obvius_fit_path{1} = 'rf-7-gf/stas-007/';
obvius_fit_path{2} = 'rf-10-gf/stas-010/';


data_paths{1} = '2007-03-27-1/data011-gdf/data011';
data_paths{2} = '2007-03-27-1/data008-gdf/data008';

base_fit_path = '/snle/lab/Experiments/Array/Analysis/2007-03-27-1/';
obvius_fit_path{1} = 'rf-11-gf/stas-011/';
obvius_fit_path{2} = 'rf-8-gf/stas-008/';

%%
data_paths{1} = '2007-03-02-1/data001-map-gdf/data001-map'; % ndf 5.1  % Very Good: good for parasols and on midgets
data_paths{2} = '2007-03-02-1/data016-gdf/data016/data016'; % ndf 0.6

base_fit_path = '/snle/lab/Experiments/Array/Analysis/2007-03-02-1/';
obvius_fit_path{1} = 'rf-1-map-gf/stas-001/';
obvius_fit_path{2} = 'rf-16-gf/stas-016/';


  
%%
data_paths{1} = '2009-04-13-1/data011/data011'; % ndf 4.0  % Very Good: good for parasols and on midgets
data_paths{2} = '2009-04-13-1/data010/data010'; % ndf 0.0

base_fit_path = '/snle/lab/Experiments/Array/Analysis/2009-04-13-1/';
obvius_fit_path{1} = 'rf-11-gf/stas-011/';
obvius_fit_path{2} = 'rf-10-gf/stas-010/';

%%


%%
num_sets = length(data_paths);
cell_types = {1,2,3,4};

for dset = 1:num_sets
    datarun{dset} = load_data(data_paths{dset});
    datarun{dset} = load_sta(datarun{dset}, 'load_sta', []);
    datarun{dset} = load_params(datarun{dset});
    datarun{dset} = load_neurons(datarun{dset});

    datarun{dset} = load_index(datarun{dset}, 'cell_type_overwrite', true);
    datarun{dset}.names.obvius_fit_path = [base_fit_path, obvius_fit_path{dset}];
    datarun{dset} = load_obvius_sta_fits(datarun{dset});
    datarun{dset} = get_sta_fits_from_obvius(datarun{dset}, cell_types);

    marks_params.thresh = 4.0;  % set the threshold for significant pixels (units sd)
    datarun{dset} = get_sta_summaries(datarun{dset}, cell_types, 'marks_params', marks_params);
    
    %datarun{dset} = compute_monitor_to_array_transformation(datarun{dset});
    datarun{dset} = load_monitor_alignment(datarun{dset});
end


types_to_plot = {1,2,3};

% generate profile plots for a chosen cell type across light levels
profiles = interdig_compute_and_plot_neighbor_rf(datarun{1},types_to_plot, 'foa', 1,...
                        'extension_factor', 3.0, 'profile_points', 125, 'distance_cutoff_factor', 1.5);

profiles = interdig_compute_and_plot_neighbor_rf(datarun{2},types_to_plot, 'foa', 2,...
                        'extension_factor', 3.0, 'profile_points', 125, 'distance_cutoff_factor', 1.5);


print(1, '~/Desktop/scotopic_overlap.pdf', '-dpdf')
print(2, '~/Desktop/photopic_overlap.pdf', '-dpdf')



% plot arbitrary profiles ontop of one another
figure(5); clf;
%plot(profiles.on_midget_x{1}, profiles.on_midget_y{1}, 'c')
hold on
%plot(profiles.on_midget_x{2}, profiles.on_midget_y{2}, 'c')
plot(profiles_fine.on_midget_x{1}, profiles_fine.on_midget_y{1}, 'b')
plot(profiles_fine.on_midget_x{2}, profiles_fine.on_midget_y{2}, 'B')

plot(profiles.on_parasol_x{1}, profiles.on_parasol_y{1}, 'r')
plot(profiles.on_parasol_x{2}, profiles.on_parasol_y{2}, 'r')

for tt = 1:length(types_to_plot)
    plot_rf_summaries(datarun{1}, {types_to_plot{tt}}, 'array', true, 'plot_fits', true, 'foa', tt)
end
print(3, '~/Desktop/on-midget.eps','-deps')


plot_rf_summaries(datarun{1}, {2}, 'array', true, 'plot_fits', true, 'foa', 1)

plot_rf_summaries(datarun{1}, {4}, 'array', true, 'plot_fits', true, 'foa', 2)
print(1, '~/Desktop/parasol.eps', '-deps')
print(2, '~/Desktop/midget.eps', '-deps')


%%
% make blurring example:



% make large RF
large_rf_params.center_radius = 10;
large_rf_params.dim = 2;
large_rf_params.x_size = 81;
large_rf_params.y_size = 81;
large_rf_params.center = [40,40];
large_rf = make_gaussian(large_rf_params);

figure(1)
imagesc(matrix_scaled_up(large_rf, 10))
colormap gray
axis equal
axis off
print(1, '~/Desktop/large_rf.pdf', '-dpdf')

% make small RF
small_rf_params.center_radius = 5;
small_rf_params.dim = 2;
small_rf_params.x_size = 81;
small_rf_params.y_size = 81;
small_rf_params.center = [40,40];
small_rf = make_gaussian(small_rf_params);

figure(2)
imagesc(matrix_scaled_up(small_rf,10))
colormap gray
axis equal
axis off
print(2, '~/Desktop/small_rf.pdf', '-dpdf')

% make blur RF
blur_rf_params.center_radius = 6;
blur_rf_params.dim = 2;
blur_rf_params.x_size = 81;
blur_rf_params.y_size = 81;
blur_rf_params.center = [40,40];
blur_rf = make_gaussian(blur_rf_params);

figure(3)
imagesc(matrix_scaled_up(blur_rf,10))
colormap gray
axis equal
axis off
print(3, '~/Desktop/blur_rf.pdf', '-dpdf')

large_scotopic_rf = filter2(blur_rf, large_rf);
small_scotopic_rf = filter2(blur_rf, small_rf);

figure(4)
imagesc(matrix_scaled_up(large_scotopic_rf,10))
colormap gray
axis equal
axis off
print(4, '~/Desktop/large_scotopic_rf.pdf', '-dpdf')

figure(5)
imagesc(matrix_scaled_up(small_scotopic_rf,10))
colormap gray
axis equal
axis off
print(5, '~/Desktop/small_scotopic_rf.pdf', '-dpdf')


%%

%%% figure out whether there is a single convolution that can explain RF size change in ON cells.

% make blur RF
parasol_rf_params.center_radius = 74.6;
parasol_rf_params.dim = 2;
parasol_rf_params.x_size = 501;
parasol_rf_params.y_size = 501;
parasol_rf_params.center = [250,250];
parasol_rf = make_gaussian(parasol_rf_params);

figure(3)
imagesc(parasol_rf)
colormap gray
axis equal
axis off

% estimage width of Gaussian
[parasol_profile_x, parasol_profile_y] = rf_profile(parasol_rf, parasol_rf_params.center);

% make blur RF
blur_rf_params.center_radius = 8.4;
blur_rf_params.dim = 2;
blur_rf_params.x_size = 401;
blur_rf_params.y_size = 401;
blur_rf_params.center = [200, 200];
blur_rf = make_gaussian(blur_rf_params);

figure(3)
imagesc(blur_rf)
colormap gray; axis equal; axis off;


% blur the parasol RF
test_rf = filter2(blur_rf, parasol_rf);
figure(3)
imagesc(test_rf)
colormap gray; axis equal; axis off;

% get the profile
[test_profile_x, test_profile_y] = rf_profile(test_rf, parasol_rf_params.center);

% fit the profile to get the new sigma
fit_params.center_radius = 90;
fit_params.center_scale = 1;
fit_params.surround_radius = 0;
fit_params.surround_scale = 0;
fit_params.fit_surround_radius = false;
fit_params.fit_surround_scale = false;
fit_profile(test_profile_x, test_profile_y, fit_params)



% now see if the same blur function will work on the midget cells

% make blur RF
midget_rf_params.center_radius = 32.6;
midget_rf_params.dim = 2;
midget_rf_params.x_size = 501;
midget_rf_params.y_size = 501;
midget_rf_params.center = [250,250];
midget_rf = make_gaussian(midget_rf_params);

imagesc(midget_rf)
colormap gray

% filter the midget RF
test_rf = filter2(blur_rf,midget_rf);

% get the profile
[test_profile_x, test_profile_y] = rf_profile(test_rf, parasol_rf_params.center);

% fit the profile to get the new sigma
fit_params.center_radius = 20;
fit_params.center_scale = 1;
fit_params.surround_radius = 0;
fit_params.surround_scale = 0;
fit_params.fit_surround_radius = false;
fit_params.fit_surround_scale = false;
fit_profile(test_profile_x, test_profile_y, fit_params)





