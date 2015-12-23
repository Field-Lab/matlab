%-- Unknown date --%
spike_times = spike_find_spikes(dataset_saved,[],45);length(spike_times)*900
help linkaxes
spike_times = spike_find_spikes(dataset_saved,[],45);length(spike_times)*900
length(spike_times)
spike_times = spike_find_spikes(dataset_saved,[],45);length(spike_times)
size(filter2(visionFilter,data(ee,:),'same'))
ee
size(filter2(visionFilter,data(2,:),'same'))
size(filter2(visionFilter,data(3,:),'same'))
size(filter2(visionFilter,data(4,:),'same'))
size(filter2(visionFilter,data(5,:),'same'))
size(filter2(visionFilter,data(6,:),'same'))
size(filter2(visionFilter,data(7,:),'same'))
num_samples
chunk_end
chunk_begin
size(data)
size(rawFile.getData(dataset.electrodes(ee),1000,5))
chunk_end - chunk_begin
size(rawFile.getData(dataset.electrodes(ee),1000,1))
size(fdata)
size(rawFile.getData(dataset.electrodes(ee),1000,1))
whos
fdata(ee,:) = filter2(visionFilter,data(ee,:),'same');
spike_times = spike_find_spikes(dataset_saved,[],45);length(spike_times)
size(data)
chunk_end
chunk_begin
size(rawFile.getData(dataset.electrodes(ee),1,chunk_end - chunk_begin))
size(rawFile.getData(dataset.electrodes(ee),1,5))
size(rawFile.getData(dataset.electrodes(ee),0,5))
spike_times = spike_find_spikes(dataset_saved,[],45);spike_times
spike_times = spike_find_spikes(dataset_saved,[],45);spike_times'
spike_times = spike_find_spikes(dataset_saved,[],45);length(spike_times)
spike_times'
spike_times = spike_find_spikes(dataset_saved,[],45);length(spike_times)
size(rawFile.getData(dataset.electrodes(ee),chunk_begin + 1,chunk_end - chunk_begin))
chunk_begin
size(rawFile.getData(dataset.electrodes(ee),chunk_begin,chunk_end - chunk_begin))
spike_times = spike_find_spikes(dataset_saved,[],45);length(spike_times)
size(rawFile.getData(dataset.electrodes(ee),chunk_begin,chunk_end - chunk_begin))
spike_times = spike_find_spikes(dataset_saved,[],45);length(spike_times)
spike_times'
spike_times = spike_find_spikes(dataset_saved,[],45);length(spike_times)
%-- 7/31/07  6:14 PM --%
load('/Volumes/War/Data/Machado/matlab/saved datasets/2005-04-26-1-data006-93.mat'); new_spike_times = spike_find_spikes(dataset_saved,[],99);
%-- 7/31/07  6:33 PM --%
spike_cluster(1, 2)
load('/Volumes/Brokedown/Data/Gauthier/matlab/saved datasets/2005-04-26-0-data002-276.mat')
whos
dataset_saved
load('/Volumes/Brokedown/Data/Gauthier/matlab/saved projections/2005-04-26-0-data002-276-pca.mat')
load('/Volumes/Brokedown/Data/Gauthier/matlab/saved projections/2005-04-26-0-data002-276-lda.mat')
load('/Volumes/Brokedown/Data/Gauthier/matlab/saved projections/2005-04-26-0-data002-276-nwpca-center.mat')
clf
clc
load('/Volumes/Brokedown/Data/Gauthier/matlab/saved datasets/2005-04-26-0-data002-99.mat')
new_spike_times = spike_find_spikes(dataset,params,figure_to_use)
new_spike_times = spike_find_spikes(dataset_saved,[],13)
who
dataset_saved
load('/Volumes/War/Data/Machado/matlab/saved datasets/2005-04-26-1-data006-99.mat');
dataset_saved
%-- 8/1/07  1:35 PM --%
a.sa = 1 2 3
a.sa = [1 2 3]
a.b = [11 11]
single(a)
single(sa.a)
single(a.sa)
155 + ans
whos
a
a.sa
single(a.sa) + double(a.sa)
load('/Volumes/War/Data/Machado/matlab/saved datasets/2005-04-26-1-data006-99.mat');
new_spike_times = spike_find_spikes(dataset_saved,[],99);
[dataset] = spike_load_dataset(dataset, spike_times);
[dataset_saved] = spike_load_dataset(dataset_saved, new_spike_times);
[dataset_saved] = spike_load_dataset(dataset_saved, new_spike_times');
clear all
load('/Volumes/War/Data/Machado/matlab/saved datasets/2005-04-26-1-data006-99.mat');
load('/private/var/automount/netapp/snle/home/tamachado/Desktop/spiketimes.mat')
[dataset_saved] = spike_load_dataset(dataset_saved, new_spike_times');
dataset_saved
dataset_saved.raw_spikes = single(dataset_saved.raw_spikes)
[dataset_saved] = spike_load_dataset(dataset_saved, new_spike_times');
whos
new_spike_times = single(new_spike_times)
new_spike_times = single(new_spike_times);
[dataset_saved] = spike_load_dataset(dataset_saved, new_spike_times');
dataset_saved = spike_load_dataset(dataset_saved, new_spike_times');
dataset_saved = spike_load_dataset(dataset_saved, new_spike_times(1:100000)');
dataset_saved = spike_plot_and_select_spikes(dataset_saved, 'pca' , [1:7] , 67, 2);
1
dataset_saved = spike_plot_and_select_spikes(dataset_saved, 'pca' , [1:7] , 67, 2);
1
dataset_saved = spike_plot_and_select_spikes(dataset_saved, 'pca' , [1:7] , 67, 2);
2
dataset_saved2 = spike_load_dataset(dataset_saved, new_spike_times(100000:200000)');
dataset_saved = spike_plot_and_select_spikes(dataset_saved2, 'pca' , [1:7] , 32, 2);
1
dataset_saved = spike_select_several_clusters(dataset_saved.pca{2},11,5);
1
% dataset_name = '2005-04-26-1-data006';
% num_electrodes = 512;
% raw_data_file = '/Volumes/War/Data/Machado/2005-04-26-1/data006/';
% neurons_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.neurons';
% spikes_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.spikes';
% model_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.model';
% prj_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.prj';
% mdf_file = '/snle/acquisition/mdf/RGB-10-8-0.48-11111.mdf';
% saved_datasets_location = '/Volumes/War/Data/Machado/matlab/saved datasets';
% saved_projections_location = '/Volumes/War/Data/Machado/matlab/saved projections';
% neurons_file_output_location = '/Volumes/War/Data/Machado/matlab/saved neurons';
% dataset_name = '2005-04-26-1-data006';
% num_electrodes = 512;
% raw_data_file = '/Volumes/War/Data/Machado/2005-04-26-1/data006/';
% neurons_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.neurons';
% spikes_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.spikes';
% model_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.model';
% prj_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.prj';
% mdf_file = '/snle/acquisition/mdf/RGB-10-8-0.48-11111.mdf';
% saved_datasets_location = '/Volumes/War/Data/Machado/matlab/saved datasets';
% saved_projections_location = '/Volumes/War/Data/Machado/matlab/saved projections';
% neurons_file_output_location = '/Volumes/War/Data/Machado/matlab/saved neurons';
% dataset_name = '2005-04-26-1-data006';
% num_electrodes = 512;
% raw_data_file = '/Volumes/War/Data/Machado/2005-04-26-1/data006/';
% neurons_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.neurons';
% spikes_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.spikes';
% model_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.model';
% prj_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.prj';
% mdf_file = '/snle/acquisition/mdf/RGB-10-8-0.48-11111.mdf';
% saved_datasets_location = '/Volumes/War/Data/Machado/matlab/saved datasets';
% saved_projections_location = '/Volumes/War/Data/Machado/matlab/saved projections';
% neurons_file_output_location = '/Volumes/War/Data/Machado/matlab/saved neurons';
%-- 8/2/07  5:13 PM --%
clc
load('/Volumes/War/Data/Machado/matlab/saved projections/2005-04-26-1-data006-57-pca.mat')
[nClusters, assignments, M, R, likelihood] = spike_mcmc_clustering(spike_projections, nPoints, iterations, mean_record, covariance_record, K_record)
proj_struct
[nClusters, assignments, M, R, likelihood] = spike_mcmc_clustering(proj_struct.spike_projections, 1, 1, mean_record, covariance_record, K_record)
%cd /var/automount/netapp/snle/home/gauthier/Desktop/PCA testing/final version/crp/crp;
%[class_id, mean_record, covariance_record, K_record, lP_record, alpha_record] = sampler(test_data', iterations);
%cd cwd;
%-- 8/6/07 10:53 AM --%
load('/Volumes/War/Data/Machado/matlab/saved projections/2005-04-26-1-data006-57-pca.mat')
%-- 8/7/07 10:59 AM --%
load('/Volumes/Sunshine/Data/Machado/matlab/saved projections/2005-04-26-0-data002-58-pca.mat');
proj-struct
proj_struct
[class_id, mean_record, covariance_record, K_record, lP_record, alpha_record] = sampler(proj_struct.spike_projections(100000:102000,1:2)', 1000);
[nClusters, assignments, M, R, likelihood] = spike_mcmc_clustering(proj_struct.spike_projections, 1, 1, mean_record, covariance_record, K_record);
load('/Volumes/Sunshine/Data/Machado/matlab/saved datasets/2005-04-26-1-data006-58.mat'); [clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 58);
load('/Volumes/Sunshine/Data/Machado/matlab/saved datasets/2005-04-26-0-data006-58.mat'); [clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 58);
load('/Volumes/Sunshine/Data/Machado/matlab/saved datasets/2005-04-26-0-data002-58.mat'); [clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 58);
[clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 58);
proj_struct
dataset_saved
[clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 58);
115636*2
load('/Volumes/Sunshine/Data/Machado/matlab/saved datasets/2005-04-26-0-data002-15.mat'); [clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 15);
load('/Volumes/Sunshine/Data/Machado/matlab/saved projections/2005-04-26-0-data002-15-pca.mat'); [class_id, mean_record, covariance_record, K_record, lP_record, alpha_record] = sampler(proj_struct.spike_projections(100000:102000,1:2)', 1000);
load('/Volumes/Sunshine/Data/Machado/matlab/saved projections/2005-04-26-0-data002-15-pca.mat'); [class_id, mean_record, covariance_record, K_record, lP_record, alpha_record] = sampler(proj_struct.spike_projections(100000:102000,1:2)', 500);
%-- 8/7/07  2:54 PM --%
load('/Volumes/Sunshine/Data/Machado/matlab/saved projections/2005-04-26-0-data002-58-pca.mat');
[class_id, mean_record, covariance_record, K_record, lP_record, alpha_record] = sampler(proj_struct.spike_projections(1:2000,1:2)', 500);
load('/Volumes/Sunshine/Data/Machado/matlab/saved datasets/2005-04-26-1-data006-58.mat'); [clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 58);
load('/Volumes/Sunshine/Data/Machado/matlab/saved datasets/2005-04-26-0-data002-58.mat'); [clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 58);
dataset_saved
size(dataset_saved.raw_spikes,1)
[clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 58);
load('/Volumes/Sunshine/Data/Machado/matlab/saved projections/2005-04-26-0-data002-9-nwpca-center.mat');
[class_id, mean_record, covariance_record, K_record, lP_record, alpha_record] = sampler(proj_struct.spike_projections(1:2000,1:2)', 500);
load('/Volumes/Sunshine/Data/Machado/matlab/saved datasets/2005-04-26-0-data002-9.mat'); [clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 9);
[clusters, sorter, assignments] = spike_mcmc_clustering(proj_struct, dataset_saved, 1, 1000, mean_record, covariance_record, K_record, 9);
%-- 8/8/07 12:59 PM --%
%-- 8/13/07  2:20 PM --%
bitshift(10,2)
bitshift(2,2)
bitshift(1,2)
%-- 8/17/07  5:54 PM --%
clear java
whos
model
model.getSpikeTimes()
clear java
model.getSpikeTimes()
model.loadElectrode(1)
data.getSpikeTimes
whos
blah = neuronsData(prj_file)
blah = neuronsData(sprintf(prj_file))
blah = edu.salk.snl.cellfinder.data.NeuronsData(sprintf(prj_file))
blah.getSpikeTimes()
clusterData = blah.loadElectrode(1)
clusterData.getSpikeTimes()
spikeTimes  = clusterData.getSpikeTimes()
spikeTimes  = clusterData.getSpikeTimes();
size(spikeTimes)
st = spikeTimes(spikeTimes > 0);
st
a = []
%-- 9/24/07  9:58 AM --%
1105554/106598
d = ans;
whos
spike_times
neuronsData = edu.salk.snl.cellfinder.data.NeuronsData(dataset.projections_file)
neuronsData = edu.salk.snl.cellfinder.data.NeuronsData('/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.prj')
neuronsData.loadElectrode(1);
st = cd.getSpikeTimes();
cd = neuronsData.loadElectrode(1);
st = cd.getSpikeTimes();
st
st = st > 0;
st
st = cd.getSpikeTimes();
st = st(st > 0);
st
size(st)
figure(99); plot(st);
figure(15); hold on; plot(st, 'g');
cd = neuronsData.loadElectrode(2);
st = cd.getSpikeTimes();
st = st(st > 0);
figure(15); hold on; plot(st, 'g');
cd = neuronsData.loadElectrode(3);
st = cd.getSpikeTimes();
st = st(st > 0);
figure(15); hold on; plot(st, 'g');
cd = neuronsData.loadElectrode(4);
st = cd.getSpikeTimes();
st = st(st > 0);
figure(15); hold on; plot(st, 'g');
cd = neuronsData.loadElectrode(5);
st = cd.getSpikeTimes();
st = st(st > 0);
figure(15); hold on; plot(st, 'g');
%-- 9/28/07 10:53 AM --%
dataset_name = '2005-09-13-1-data000';
num_electrodes = 512;
raw_data_file = '/Volumes/Creampuff/Data/gfield/2005-09-13-1/data000/';
neurons_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.neurons';
spikes_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.spikes';
model_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.model';
prj_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.prj';
mdf_file = '/snle/acquisition/mdf/RGB-10-8-0.48-11111.mdf';
saved_datasets_location = '/Volumes/Creampuff/Data/Machado/matlab/saved datasets';
saved_projections_location = '/Volumes/Creampuff/Data/Machado/matlab/saved projections';
neurons_file_output_location = '/Volumes/Creampuff/Data/Machado/matlab/saved neurons';
dataset_name = '2005-09-13-1-data000';
num_electrodes = 512;
raw_data_file = '/Volumes/Creampuff/Data/gfield/2005-09-13-1/data000/';
neurons_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.neurons';
spikes_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.spikes';
model_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.model';
prj_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.prj';
mdf_file = '/snle/acquisition/mdf/RGB-10-8-0.48-11111.mdf';
saved_datasets_location = '/Volumes/Creampuff/Data/Machado/matlab/saved datasets';
saved_projections_location = '/Volumes/Creampuff/Data/Machado/matlab/saved projections';
neurons_file_output_location = '/Volumes/Creampuff/Data/Machado/matlab/saved neurons';
dataset_name = '2005-09-13-1-data000';
num_electrodes = 512;
raw_data_file = '/Volumes/Creampuff/Data/gfield/2005-09-13-1/data000/';
neurons_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.neurons';
spikes_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.spikes';
model_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.model';
prj_file = '/Volumes/Creampuff/Analysis/Machado/2005-09-13-1/data000/data000.prj';
mdf_file = '/snle/acquisition/mdf/RGB-10-8-0.48-11111.mdf';
saved_datasets_location = '/Volumes/Creampuff/Data/Machado/matlab/saved datasets';
saved_projections_location = '/Volumes/Creampuff/Data/Machado/matlab/saved projections';
neurons_file_output_location = '/Volumes/Creampuff/Data/Machado/matlab/saved neurons';
clc
load('/Data/Machado/matlab/saved datasets/2005-09-27-4-data000-2.mat')
dd = 1;dataset{dd} = dataset_saved;clear dataset_saved;
dataset{dd}
dataset{dd} = spike_plot_and_select_spikes(dataset{dd},'pca',[1:7],[3],1);
1
load('/Data/Machado/matlab/saved datasets/2005-09-27-4-data000-295.mat')
dd = 1;dataset{dd} = dataset_saved;clear dataset_saved;
dataset{dd} = spike_plot_and_select_spikes(dataset{dd},'pca',[1:7],[3],1);
1
dataset{dd}.pca{1}
spikes = dataset{dd}.pca{1}.spike_projections;
load('/Volumes/War/Data/Machado/matlab/saved datasets/2007-04-27-data003-1.mat')
dd = 1;dataset{dd} = dataset_saved;clear dataset_saved;
dataset{dd} = spike_plot_and_select_spikes(dataset{dd},'pca',[1:7],[3],1);
%-- 9/28/07  1:29 PM --%
load('/Volumes/War/Data/Machado/matlab/saved datasets/2007-04-27-data003-1.mat')
dd = 1;dataset{dd} = dataset_saved;clear dataset_saved;
dataset{dd} = spike_plot_and_select_spikes(dataset{dd},'pca',[1:7],[3],1);
1
state = 634;
rand('state', state);
%get random points
randNums = unidrnd(nSpikes, 1, nPoints);
%assign the selected random points
test_data = zeros(nPoints, length(DIMS));
for point = 1:nPoints
test_data(point,:) = spike_projections(randNums(point),:);
end
nSpikes = 10000
nPoints = 1000
state = 634;
rand('state', state);
%get random points
randNums = unidrnd(nSpikes, 1, nPoints);
%assign the selected random points
test_data = zeros(nPoints, length(DIMS));
for point = 1:nPoints
test_data(point,:) = spike_projections(randNums(point),:);
end
DIMS = 1:5
state = 634;
rand('state', state);
%get random points
randNums = unidrnd(nSpikes, 1, nPoints);
%assign the selected random points
test_data = zeros(nPoints, length(DIMS));
for point = 1:nPoints
test_data(point,:) = spike_projections(randNums(point),:);
end
spikes = dataset{dd}.pca{1}.spike_projections;
spike_projections = spikes
spike_projections = spikes;
state = 634;
rand('state', state);
%get random points
randNums = unidrnd(nSpikes, 1, nPoints);
%assign the selected random points
test_data = zeros(nPoints, length(DIMS));
for point = 1:nPoints
test_data(point,:) = spike_projections(randNums(point),:);
end
randNums
state = 634;
rand('state', state);
%get random points
randNums = unidrnd(nSpikes, 1, nPoints);
%assign the selected random points
test_data = zeros(nPoints, length(DIMS));
for point = 1:nPoints
test_data(point,:) = spike_projections(randNums(point),:);
end
%-- 10/1/07 11:56 AM --%
%-- 10/2/07  4:02 PM --%
proj_struct
proj_struct.dimensions_to_plot
proj_struct
clear java
java
javaclasspath
clear java
proj_struct
mean(spike_projections(:,1:5))
mean(proj_struct.spike_projections(:,1:5))
cov(proj_struct.spike_projections(:,1:5))
3 * cov(proj_struct.spike_projections(:,1:5))
mean_record
3 * cov(proj_struct.spike_projections(:,1:5))
mean(proj_struct.spike_projections(:,1:5))
[1 2 3 4]
a = ans;
a
b = [1 2 3 4; 2 3 4 5];
[a; b]
proj_struct
%-- 10/3/07  6:11 PM --%
dataset_name = '2005-04-26-1-data006';
num_electrodes = 512;
raw_data_file = '/Volumes/War/Data/Machado/2005-04-26-1/data006/';
neurons_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.neurons';
spikes_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.spikes';
model_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.model';
prj_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.prj';
mdf_file = '/snle/acquisition/mdf/RGB-10-8-0.48-11111.mdf';
saved_datasets_location = '/Volumes/War/Data/Machado/matlab/saved datasets';
saved_projections_location = '/Volumes/War/Data/Machado/matlab/saved projections';
neurons_file_output_location = '/Volumes/War/Data/Machado/matlab/saved neurons';
dataset_name = '2005-04-26-1-data006';
num_electrodes = 512;
raw_data_file = '/Volumes/War/Data/Machado/2005-04-26-1/data006/';
neurons_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.neurons';
spikes_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.spikes';
model_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.model';
prj_file = '/Volumes/War/Analysis/Machado/2005-04-26-1/data006/data006.prj';
mdf_file = '/snle/acquisition/mdf/RGB-10-8-0.48-11111.mdf';
saved_datasets_location = '/Volumes/War/Data/Machado/matlab/saved datasets';
saved_projections_location = '/Volumes/War/Data/Machado/matlab/saved projections';
neurons_file_output_location = '/Volumes/War/Data/Machado/matlab/saved neurons';
disp(lasterr)
help lasterr
disp(lasterror)
lasterror
e = lasterror
e.stack
disp(e.stack)
e
e = lasterr
e = lasterror
e{stack}
e
e.stack
s = e.stack
size(s)
s(1)
disp(s)
e.stack
%-- 10/5/07  3:12 PM --%
%-- 10/8/07  2:30 PM --%
%-- 10/11/07  4:27 PM --%
