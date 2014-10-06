
if ~exist('rawFile','var')
    rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile('/Data/Gauthier/2007-09-18-4/data002');
    datarun = load_data('2007-09-18-4/data002-nwpca/data002');
    datarun.piece.array_id = 500;
    datarun = load_neurons(datarun,'sort_cell_ids',true);
    datarun = load_ei(datarun,datarun.cell_ids(1));
end




% find closest spikes
if 0
    % load ISI(t) and ISI(t-1)
    xx=spikes(3:end)-spikes(2:end-1);yy=spikes(2:end-1)-spikes(1:end-2);
    % plot
    figure;plot(xx,yy,'.')
    % find three spikes very close together
    disp(spikes(xx<500 & yy<1000))
end



cell_id = 1326;
electrodes = [89 296 483];
plot_start = 53127629-200;
plot_duration = 2000;
y_lim = [-300 300];
multiplier = [1 20 5];


% get data
samplingRate = edu.ucsc.neurobiology.vision.util.VisionParams.samplesPerMillisecond * 1000;
data = rawFile.getData(plot_start, plot_duration);
disp(size(data))

% get all spike times
spikes = datarun.spikes{get_cell_indices(datarun,cell_id)} * samplingRate;

% identify spike times relevant to this region
spikes_here = spikes((spikes>plot_start) & (spikes<(plot_start+plot_duration))) - plot_start;
disp(size(spikes_here))

% get ei
ei = get_ei(datarun,cell_id);

% plot
figure(1);clf
for ee=1:length(electrodes)
    % plot raw voltage
    subplot(length(electrodes),2,2*(ee-1)+1);hold on
    plot(1:plot_duration, data(:,electrodes(ee)+1) - mean(data(:,electrodes(ee)+1)))
    for ss=1:length(spikes_here)
        plot([1 1]*spikes_here(ss),y_lim,'r')
    end
    ylim(y_lim)
    
    % plot ei
    subplot(length(electrodes),2,2*ee)
    plot(ei(electrodes(ee),:)*multiplier(ee))
    title(num2str(multiplier(ee)))
    ylim(y_lim)
end

%linkaxes(get(1,'child'))