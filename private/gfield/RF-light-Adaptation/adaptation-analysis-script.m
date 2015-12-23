% RF analysis for Adaptation project

%%
data_paths{1} = '2007-02-06-4/data001-mg/data001/data001'; % ndf 5.1
data_paths{2} = '2007-02-06-4/data004-mg/data004/data004'; % ndf 3.9
data_paths{3} = '2007-02-06-4/data005-mg/data005/data005'; % ndf 3.0
data_paths{4} = '2007-02-06-4/data007-mg/data007/data007'; % ndf 2.3
data_paths{5} = '2007-02-06-4/data009-mg/data009/data009'; % ndf 2.0
data_paths{6} = '2007-02-06-4/data011-mg/data011/data011'; % ndf 1.6
data_paths{7} = '2007-02-06-4/data013-mg/data013/data013'; % ndf 1.3
data_paths{8} = '2007-02-06-4/data015-mg/data015/data015'; % ndf 0.6

base_fit_path = '/snle/lab/Experiments/Array/Analysis/2007-02-06-4/';

obvius_fit_path{1} = 'rf-1-gf/stas-001/';
obvius_fit_path{2} = 'rf-4-gf/stas-004/';
obvius_fit_path{3} = 'rf-5-gf/stas-005/';
obvius_fit_path{4} = 'rf-7-gf/stas-007/';
obvius_fit_path{5} = 'rf-9-gf/stas-009/';
obvius_fit_path{6} = 'rf-11-gf/stas-011/';
obvius_fit_path{7} = 'rf-13-gf/stas-013/';
obvius_fit_path{8} = 'rf-15-gf/stas-015/';

num_sets = length(data_paths);
cell_types = {1,2,3,4};
%%

%% 
data_paths{1} = '2007-03-08-4/data003/data003'; % ndf 2.0
data_paths{1} = '2007-03-08-4/data005/data005'; % ndf 0.0

%%


for dset = 1:num_sets
    datarun{dset} = load_data(data_paths{dset});
    datarun{dset} = load_sta(datarun{dset}, 'load_sta', []);
    datarun{dset} = load_params(datarun{dset});
    datarun{dset} = load_neurons(datarun{dset});
    datarun{dset} = load_index(datarun{dset}, 'cell_type_overwrite', true);
    datarun{dset}.names.obvius_fit_path = [base_fit_path, obvius_fit_path{dset}];
    datarun{dset} = load_obvius_sta_fits(datarun{dset});
    datarun{dset} = compute_monitor_to_array_transformation(datarun{dset});

    datarun{dset} = get_sta_fits_from_obvius(datarun{dset}, cell_types);

    marks_params.thresh = 4.0;  % set the threshold for significant pixels (units sd)
    datarun{dset} = get_sta_summaries(datarun{dset}, cell_types, 'marks_params', marks_params);
end

plot_rf_summaries(datarun{dset}, {1}, 'foa', dset, 'plot_fits', true, 'array', true)


% calculate the RF radii for the specified cell type

cell_type = {1};
on_parasol_mean_radius = zeros(num_sets,1);
for dset = 1:num_sets
    temp_radii = get_rf_fit_radius(datarun{dset}, cell_type,...
                    'fits_to_use', 'obvius', 'units', 'microns', 'microns_per_pixel', 5.5);
    on_parasol_mean_radius(dset) = mean(temp_radii);
    plot_rf_summaries(datarun{dset}, cell_type, 'foa', dset, 'plot_fits', true, 'array', true)

end

cell_type = {2};
off_parasol_mean_radius = zeros(num_sets,1);
for dset = 1:num_sets
    temp_radii = get_rf_fit_radius(datarun{dset}, cell_type,...
                    'fits_to_use', 'obvius', 'units', 'microns', 'microns_per_pixel', 5.5);
    plot_rf_summaries(datarun{dset}, cell_type, 'foa', dset, 'plot_fits', true, 'array', true)
end

cell_type = {3};
on_midget_mean_radius = zeros(num_sets,1);
for dset = 1:num_sets
    temp_radii = get_rf_fit_radius(datarun{dset}, cell_type,...
                    'fits_to_use', 'obvius', 'units', 'microns', 'microns_per_pixel', 5.5);
    plot_rf_summaries(datarun{dset}, cell_type, 'foa', dset, 'plot_fits', true, 'array', true)
end


cell_type = {4};
off_midget_mean_radius = zeros(num_sets,1);
for dset = 1:num_sets
    temp_radii = get_rf_fit_radius(datarun{dset}, cell_type,...
                    'fits_to_use', 'obvius', 'units', 'microns', 'microns_per_pixel', 5.5);
    off_midget_mean_radius(dset) = mean(temp_radii);
    plot_rf_summaries(datarun{dset}, cell_type, 'foa', dset, 'plot_fits', true, 'array', true)
end


figure(1); clf;
plot(on_parasol_mean_radius, 'ko-')
hold on
plot(off_parasol_mean_radius, 'ro-')
plot(on_midget_mean_radius,'ks-')
plot(off_midget_mean_radius,'rs-')
hold off


% generate profile plots for a chosen cell type across light levels
for dset =2:num_sets
    profiles = interdig_compute_and_plot_neighbor_rf(datarun{dset},{1,2,3,4}, 'foa', dset);
end








