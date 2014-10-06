clear all

movieNo = 4;
%xLimits = [0 1];
xLimits = [];
orderNo = 3; %use 0 to signify simultaneous stimulation and -1 to signify isolated stimulation
offset = 5; %1 2 or 5 (ms)


pathToData = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data011';


pxPerMs = 50;
rasterMarkerSize = 3;


% dataset-specific parameters

neuronIDs = [1 228 332 378]; %targets should correspond with pattern numbers

electrodes = [60 14 21 26]; %stimulating electrodes, in order corresponding with target neurons in neuronIDs

%orders in which corresponding patterns 
%(cell array index refers to which sequence is applied, vector index corresponds to 'electrodes', value = ith pulse applied in sequence) are played
patternOrders = {[4 1 3 2];
                 [1 4 3 2];
                 [1 3 2 4];
                 [1 2 3 4]};
           

spatiotempPatternsFigureMakerNew(pathToData, neuronIDs, electrodes, patternOrders, offset, orderNo, movieNo,...
    'xLimits', xLimits, 'pxPerMs', pxPerMs, 'rasterMarkerSize', rasterMarkerSize)


keyboard

%% plot mosaic and auto/cross correlograms

datarun = load_data('2011-06-24-5/data009/data009');

datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun, 'load_sta', []);
%datarun = load_ei(datarun, []);
datarun = load_neurons(datarun);
datarun = get_sta_fits_from_vision(datarun, 'all');

% mosaics plot
figure
for ii = 1:4
    plot_rf_summaries(datarun, neuronIDs(ii), 'fit_color', colorScheme(ii,:))
end

return

figure('position', [100 100 800 800])
correlation_axes = cell(length(neuronIDs), length(neuronIDs));
options = struct('offset', 0.02 ,'scale','ms','shuffle','none', 'dt', 0.01/32);

for ii = 1:4
    for jj = 1:4
        correlation_axes{ii,jj} = axes('position', [0.22*(4-jj) + 0.1, 0.22*(4-ii) + 0.1, 0.15, 0.15]);
        plot_ccf(datarun, [neuronIDs(ii) neuronIDs(jj)], 'fig_or_axes', correlation_axes{ii,jj}, 'verbose', false, 'options', options)
        title([num2str(neuronIDs(ii)) ' x ' num2str(neuronIDs(jj))])
    end
end


%% plot eis

[xCoords yCoords] = getElectrodeCoords61();
pathToEi = '/snle/lab/Experiments/Array/Analysis/2009-09-07-0/data010/data010.ei';

nContours = 6;

figure
hold on

plot(xCoords, yCoords, 'ok','MarkerSize', 2, 'MarkerFaceColor', 'k')

for i = 1:length(neuronIDs)
    neuronID = neuronIDs(i);
    % get neuron's average spike waveforms from VISION analysis
    eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
    ei = eiFile.getImage(neuronID); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples
    clear eiFile

    % calculate maximum waveform value on each electrode (absolute value)
    eiAmps = zeros(64, 1);
    for j = 1:64
        if ~(j==9||j==25||j==57)
            eiAmps(j) = max(max(abs(ei(1,j+1,:))));
        end
    end
    eiAmps = eiAmps/max(eiAmps);

    
    contours = hex_contour(xCoords, yCoords, eiAmps, nContours, 'fig_or_axes', [], 'contourSpacing', 'linear');
    %title(['neuron' num2str(neuronIDs(i))])
    for j = 1:nContours
        for k = 1:length(contours(j).paths)
            plot(contours(j).paths{k}(:,1), contours(j).paths{k}(:,2), 'color', colorScheme(i,:))
        end
    end
    
end

%plot asterisks on stim elecs
for i = 1:length(electrodes)
    plot(xCoords(electrodes(i)), yCoords(electrodes(i)), '*')
end

hold off
axis equal

