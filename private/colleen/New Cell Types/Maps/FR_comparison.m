clear
match = crossIdentifyNeuronIDs('/Volumes/Analysis/2016-02-17-6/streamed/data007','/Volumes/Analysis/2016-02-17-6/data015', [],[], 1,1 )

datarun.names.rrs_neurons_path = '/Volumes/Analysis/2016-02-17-6/streamed/data007_colleen/data007.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2016-02-17-6/streamed/data007_colleen/data007.params';

datarun = load_neurons(datarun);
datarun = load_params(datarun);

cell_indices = get_cell_indices(datarun, 'ON parasol');
cell_ids = get_cell_ids(datarun, 'ON parasol');

for i = 1:length(cell_ids)
match2(i,:) = match(match(:,1) == cell_ids(i), :);
end

FR_on(:,1) = match2(:,1);
FR_on(:,3) = match2(:,2);
for i = 1:length(cell_ids)
%     if ~isnan(match(match(:,1) == cell_ids(i),2));
            FR_on(i,2) = length(datarun.spikes{cell_indices(i)})./1800;
%     else
%                     FR_on(i,1) = nan;
%     end
end

clear datarun
datarun.names.rrs_neurons_path = '/Volumes/Analysis/2016-02-17-6//data015/data015.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2016-02-17-6/data015/data015.params';

datarun = load_neurons(datarun);
datarun = load_params(datarun);

cell_indices = get_cell_indices(datarun, 'ON parasol');
cell_ids = get_cell_ids(datarun, 'ON parasol');

for i = 1:length(match2(:,2))

    if ~isnan(match2(i,2))
                cell_indices = get_cell_indices(datarun, match2(i,2));

            FR_on(i,4) = length(datarun.spikes{cell_indices})./600;
    else
                    FR_on(i,4) = nan;
    end
end

FR_on(:,5) = (FR_on(:,2) - FR_on(:,4))./FR_on(:,2)*100;

%%


clear datarun
match = crossIdentifyNeuronIDs('/Volumes/Analysis/2016-02-17-6/streamed/data007','/Volumes/Analysis/2016-02-17-6/data016', [],[], 1,1 )

datarun.names.rrs_neurons_path = '/Volumes/Analysis/2016-02-17-6/streamed/data007_colleen/data007.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2016-02-17-6/streamed/data007_colleen/data007.params';

datarun = load_neurons(datarun);
datarun = load_params(datarun);

cell_indices = get_cell_indices(datarun, 'OFF parasol');
cell_ids = get_cell_ids(datarun, 'OFF parasol');

for i = 1:length(cell_ids)
match2(i,:) = match(match(:,1) == cell_ids(i), :);
end

FR_off(:,1) = match2(:,1);
FR_off(:,3) = match2(:,2);
for i = 1:length(cell_ids)
%     if ~isnan(match(match(:,1) == cell_ids(i),2));
            FR_off(i,2) = length(datarun.spikes{cell_indices(i)})./1800;
%     else
%                     FR_on(i,1) = nan;
%     end
end

clear datarun
datarun.names.rrs_neurons_path = '/Volumes/Analysis/2016-02-17-6//data016/data016.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2016-02-17-6/data016/data016.params';

datarun = load_neurons(datarun);
datarun = load_params(datarun);

cell_indices = get_cell_indices(datarun, 'ON parasol');
cell_ids = get_cell_ids(datarun, 'ON parasol');

for i = 1:length(match2(:,2))

    if ~isnan(match2(i,2))
                cell_indices = get_cell_indices(datarun, match2(i,2));

            FR_off(i,4) = length(datarun.spikes{cell_indices})./600;
    else
                    FR_off(i,4) = nan;
    end
end

FR_off(:,5) = (FR_off(:,2) - FR_off(:,4))./FR_off(:,2)*100;

f = figure;
hold on
[nb,xb] =hist(FR_on(:,5));
createPatches(xb,nb,[1 1 1],3,0.05);
[nb,xb] =hist(FR_off(:,5));
createPatches(xb,nb,[0 0 0],3,0.5);
legend('ON', 'OFF')
title({'2016-02-17-6'; ['ON parasol: n = ', num2str(size(FR_on(~isnan(FR_on(:,5)),1),1))];[ 'OFF parasol: n = ', num2str(size(FR_off(~isnan(FR_off(:,5)),1),1))]})
xlabel('% decrease in FR with gaussian mask stimulus centered on each parasol')
ylabel('Number of cells')

[~,git_hash_string] = system('git rev-parse HEAD');
fprintf('Git Hash: %s \n', git_hash_string);