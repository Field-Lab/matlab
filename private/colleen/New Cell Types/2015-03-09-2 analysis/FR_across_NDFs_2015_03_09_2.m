for dataset = 1:4
    clearvars datarun
    dataset_id = {'04', '05', '11', '19'};
    id = char(dataset_id{dataset});
file_name = ['2015-03-09-2/d04_05_11_15_19-norefit/data0', id, '-from-data004_data005_data011_data015_data019/data0', id, '-from-data004_data005_data011_data015_data019'];
% file_name = ['2015-03-09-2/data023-from-data019/data023-from-data019'];
% file_name = ['2015-03-09-2/data015-from-data019/data015-from-data004_data005_data011_data015_data019'];

datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];

% datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];
% mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111.xml';
% interpolate = false;
% cell_specification = [502,860,1024,1130,2076,2361,2618,2705,3022,3172,3213,3559,4022,4071,4238,4774,4852,5496,6518,6533,6860,7279,7671];
cell_type = {'on parasol'};


opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',0, 'load_sta', 0, 'load_sta_params', 0, 'load_all',0);
% opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 11:30% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

for i = 1:size(datarun.cell_types,2)
    result(i) = strcmp(datarun.cell_types{i}.name, 'ON parasol');
end
ON_parasol_ind = find(result == 1);
cell_specification = datarun.cell_types{ON_parasol_ind}.cell_ids;
matching = full(datarun.cell_nums);
cell_ids = matching(cell_specification);


for j =1:length(cell_ids)
num_spikes(j) = size(datarun.spikes{cell_ids(j)},1);
end
avg_FR_parasol(dataset) = mean(num_spikes);

end
clearvars datarun
file_name = ['2015-03-09-2/data015-from-data019/data015-from-data004_data005_data011_data015_data019'];
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];

% datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];
% mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111.xml';
% interpolate = false;
% cell_specification = [502,860,1024,1130,2076,2361,2618,2705,3022,3172,3213,3559,4022,4071,4238,4774,4852,5496,6518,6533,6860,7279,7671];
cell_type = {'on parasol'};


opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',0, 'load_sta', 0, 'load_sta_params', 0, 'load_all',0);
% opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 11:30% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

for i = 1:size(datarun.cell_types,2)
    result(i) = strcmp(datarun.cell_types{i}.name, 'ON parasol');
end
ON_parasol_ind = find(result == 1);
cell_specification = datarun.cell_types{ON_parasol_ind}.cell_ids;
matching = full(datarun.cell_nums);
cell_ids = matching(cell_specification);


for j =1:length(cell_ids)
num_spikes(j) = size(datarun.spikes{cell_ids(j)},1);
end
avg_FR_parasol(5) = mean(num_spikes);

clearvars datarun
file_name = ['2015-03-09-2/data023-from-data019/data023-from-data019'];
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];

% datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];
% mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111.xml';
% interpolate = false;
% cell_specification = [502,860,1024,1130,2076,2361,2618,2705,3022,3172,3213,3559,4022,4071,4238,4774,4852,5496,6518,6533,6860,7279,7671];
cell_type = {'on parasol'};


opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',0, 'load_sta', 0, 'load_sta_params', 0, 'load_all',0);
% opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 11:30% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

for i = 1:size(datarun.cell_types,2)
    result(i) = strcmp(datarun.cell_types{i}.name, 'ON parasol');
end
ON_parasol_ind = find(result == 1);
cell_specification = datarun.cell_types{ON_parasol_ind}.cell_ids;
matching = full(datarun.cell_nums);
cell_ids = matching(cell_specification);


for j =1:length(cell_ids)
num_spikes(j) = size(datarun.spikes{cell_ids(j)},1);
end
avg_FR_parasol(6) = mean(num_spikes);

ndf_order = [4,5,3,1,2,0]
time = [900,1800,1800,900,900,900];
avg_FR_time = avg_FR_parasol./time;
data_order = [4,5,11,15,19,23];

figure; plot(ndf_order, avg_FR_time, 'o')

avg_reshaped = [avg_FR_time(6) avg_FR_time(4) avg_FR_time(5) avg_FR_time(3) avg_FR_time(1) avg_FR_time(2)]

matrix = [2.40111111111111	4.34777777777778	3.07555555555556	2.47277777777778	6.69444444444444	2.98277777777778
3.35222222222222	1.27666666666667	4.02777777777778	4.40333333333333	8.17	0.471111111111111
2.81333333333333	4.16555555555556	2.99777777777778	4.32055555555556	8.24111111111111	1.10722222222222
0.638888888888889	0.79	1.27111111111111	1.50277777777778	3.74	0.522222222222222
4.24777777777778	1.73	3.54444444444444	3.45333333333333	5.87444444444444	0.540555555555556
5.13111111111111	1.66111111111111	4.01222222222222	4.15777777777778	6.48888888888889	0.405555555555556
0.84	3.13555555555556	3.93777777777778	2.23111111111111	9.96888888888889	2.20833333333333
1.72111111111111	0.33	1.32666666666667	1.45722222222222	3.39777777777778	0.13
5.57111111111111	0.765555555555556	3.50333333333333	3.82166666666667	5.12111111111111	0.0427777777777778
6.14111111111111	1.84888888888889	8.18111111111111	7.51611111111111	9.67333333333333	1.83222222222222
1.75222222222222	1.72	2.28777777777778	2.34111111111111	3.68111111111111	0.213888888888889
1.57	1.60777777777778	2.06444444444444	2.36833333333333	3.82444444444444	0.0911111111111111
23	14.0933333333333	20.2322222222222	18.8772222222222	31.8033333333333	7.80444444444444];

matrix_norm = matrix./repmat(avg_reshaped, 13,1);
x = [0,1,2,3,4,5];
x_plot = repmat(x, 13,1);
cc=jet(12);
figure
for i =1:12
plot(x, matrix_norm(i,:), 'o-', 'color', cc(i,:))
hold on
end
hold on
% plot(x, matrix(13,:), 'o-', 'Linewidth', 2 )
xlabel('NDF')
ylabel('Normalized Firing Rate Units')
title({'Dataset 2015-03-09-2'; 'Firing Rate of Large Cells Across Light Levels'; 'Normalized by ON Parasol Firing Rate'})

legend('3961', '2012', '4726', '5254', '1773', '4235', '7036', '2976', '1786', '5783', '7522', '7625')