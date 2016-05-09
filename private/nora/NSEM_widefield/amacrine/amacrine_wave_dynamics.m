
datarun = load_data('2015-10-06-0/data000-data015-norefit/data002-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015/data002-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015');
datarun_class = load_data('2015-10-06-0/data000-data015-norefit/data000-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015/data000-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015');
datarun = load_neurons(datarun);
datarun_class = load_params(datarun_class);
datarun_class = load_ei(datarun_class, 'Off Amacrine nc0');
cid = get_cell_indices(datarun_class, 'Off Amacrine nc0');
prepped_data = interleaved_data_prep(datarun, [3600*2 3600], 50, 'cell_spec', 'Off Amacrine nc0', 'datarun_class', datarun_class, 'testmovie_only', 1, 'stimulus_name', 'NSbrownian_3000_A_025.rawMovie');


%% EI approach
[x,y]= getElectrodeCoords512;
electrical_activity = zeros(512, 3700);
trial = 3;
for i_cell=1:length(cid)
    EI_cell = datarun_class.ei.eis{cid(i_cell)};
    for i_spike = 1:length(prepped_data.testspikes{trial,i_cell})
        spike_frame = floor(prepped_data.testspikes{trial,i_cell}(i_spike)*120);
        electrical_activity(:,(spike_frame:spike_frame+64)) = electrical_activity(:,(spike_frame:spike_frame+64))+EI_cell;
    end
end
for i=300:460
scatter(x,y,0.5*abs(electrical_activity(:,i))+1,'filled')
axis image
axis off
title(['Frame ' num2str(i)])
pause(0.01)
end

%% surfplot approach
spike = zeros(1800,length(cid));
centers = zeros(2, length(cid));
spikes_bin = zeros(1800,length(cid));
for i_cell=1:length(cid)
    spikes = floor(120*cell2mat(prepped_data.testspikes(:,i_cell)));
    centers(:,i_cell) = flip(datarun_class.vision.sta_fits{cid(i_cell)}.mean);
    for i = 1:1800
        spikes_bin(i,i_cell) = sum(spikes == i);
    end
end
spikes_bin = conv2(spikes_bin, gausswin(10), 'same');

%%
start = 360;
time = (0:40) + start;
plot(spikes_bin(time,:)+50*repmat(centers(1,:), [length(time),1]))

%% interp approach
spike = [0 0 0 0]; 
%centers = zeros(2, length(cid));
spikes_bin = zeros(1800,1);
for i_cell=1:length(cid)
    spikes = floor(120*cell2mat(prepped_data.testspikes(:,i_cell)));
    centers(:,i_cell) = flip(datarun_class.vision.sta_fits{cid(i_cell)}.mean);
    for i = 1:1800
        spikes_bin(i) = sum(spikes == i);
    end
    spike = [ spike; spikes_bin, centers(1)*ones(1800,1), centers(2)*ones(1800,1), (1:1800)']; 
end
spikes_bin = conv2(spikes_bin, gausswin(10), 'same');