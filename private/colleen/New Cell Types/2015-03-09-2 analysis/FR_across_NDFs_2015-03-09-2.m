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
figure; plot(ndf_order, avg_FR_time, 'o')

avg_reshaped = [avg_FR_time(6) avg_FR_time(4) avg_FR_time(5) avg_FR_time(3) avg_FR_time(1) avg_FR_time(2)]

matrix = [2161	3913	2768	4451	6025	5369
3017	1149	3625	7926	7353	848
2532	3749	2698	7777	7417	1993
575	711	1144	2705	3366	940
3823	1557	3190	6216	5287	973
4618	1495	3611	7484	5840	730
756	2822	3544	4016	8972	3975
1549	297	1194	2623	3058	234
5014	689	3153	6879	4609	77
5527	1664	7363	13529	8706	3298
1577	1548	2059	4214	3313	385
1413	1447	1858	4263	3442	164
20700	12684	18209	33979	28623	14048];

matrix_norm = matrix./repmat(avg_reshaped, 13,1);
x = [0,1,2,3,4,5];
x_plot = repmat(x, 13,1);
figure
for i =1:12
plot(x, matrix_norm(i,:), 'o-')
hold on
end
hold on
% plot(x, matrix(13,:), 'o-', 'Linewidth', 2 )
xlabel('NDF')
ylabel('Normalized Firing Rate Units')
title({'Dataset 2015-03-09-2'; 'Firing Rate of Large Cells Across Light Levels'; 'Normalized by ON Parasol Firing Rate'})

