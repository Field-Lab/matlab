clear
rawDataPath  = '/Volumes/Data/2016-04-21-8/data022';

vision_id = [4827 4824 4561 5011 4651 4981 4936 4839];
max_electrode = 322;


rawDataFile = edu.ucsc.neurobiology.vision.io.RawDataFile(rawDataPath);
elect322 = double(rawDataFile.getData(322, 0 ,100000));

datarun = load_data('/Volumes/Analysis/2016-04-21-8/data022-mVision/data022_regroup422_4818/data022/data022.neurons');
cell_id = get_cell_indices(datarun, vision_id);

for i = 1:length(cell_id)
spikes{:,i} = datarun.spikes{cell_id(i)};
spikes_first_five{:,i} = datarun.spikes{cell_id(i)}(datarun.spikes{cell_id(i)} < 5)*20000; % in sample bins

end

datarun.names.rrs_ei_path = '/Volumes/Analysis/2016-04-21-8/data022-mVision/data022_regroup422_4818/data022/data022.ei';

for i = 1:length(cell_id)
    datarun = load_ei(datarun, vision_id(i));
    waveforms(i,:) = datarun.ei.eis{cell_id(i)}(max_electrode,:);

end
figure; plot(waveforms')

figure; plot(elect322)
hold on

color = hsv(length(cell_id));
for i = 1:length(cell_id)
    plot(spikes_first_five{i}, 500*ones(length(spikes_first_five{i}),1), '.', 'MarkerSize', 15, 'MarkerFaceColor', color(i,:));
end




% w = conv(double(elect322), double(waveform));

