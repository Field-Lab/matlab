clear all

movieNo = 5;
%xLimits = [0 1];
xLimits = [];
iPatternNo = 2;


pathToData = '/snle/lab/Experiments/Array/Analysis/2009-09-07-0/data008';

colorScheme(1,:) = [27 117 187]/256;
colorScheme(2,:) = [190 30 45]/256;
colorScheme(3,:) = [41 180 115]/256;
colorScheme(4,:) = [0.8 0.5 0];

%these should match actual order of pulses, determined by corresponding pattern times!!!
patternNosAll = {1;
                [201 301 401 501];
                [202 302 402 502];
                [203 303 403 503];
                [204 304 404 504]};
             
%times at which corresponding patternNos are played (in ms)
patternTimesAll = {0;
                  [1 0 2 3];
                  [0 2 3 1];
                  [3 2 0 1];
                  [2 0 1 3]};


neuronIDs = [151 481 496 916]; %targets should correspond with pattern numbers

%electrodes only required for simultaneous stim
electrodes = [11 33 34 59]; %stimulating electrodes, in order corresponding with target neurons in neuronIDs

%neurons and patterns will be displayed from bottom to top (to reverse, reverse patternTimes)
spatiotempPatternsFigureMaker(pathToData, neuronIDs, patternNosAll{iPatternNo}, patternTimesAll{iPatternNo}, movieNo, 'xLimits', xLimits, 'electrodes', electrodes)

%% plot mosaic and auto/cross correlograms

datarun = load_data('2009-09-07-0/data010/data010');

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

