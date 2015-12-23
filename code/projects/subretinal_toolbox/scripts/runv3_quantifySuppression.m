clear;

% 2014-01-15-0: RCS
statisticsActivationFolder = '/media/MEA_PROCESSED_8/2014-01-15-0/data/data001/statistics/activation';
neuronIDList = [152,437,676,721,738,781,1817,4081,4351,4381,4486,4532,4923,4951,5326,5446,5521,5598,6167,7111,7201];
neuronIDType = [  2,  2,  2,  1,  1,  1,   4,   2,   2,   2,   2,   1,   1,   1,   2,   2,   2,   4,   3,   2,   2];

% 2014-01-17-0: RCS
statisticsActivationFolder = '/media/MEA_PROCESSED_8/2014-01-17-0/data/data000/statistics/activation';
neuronIDList = [3,77,211,241,436,556,587,601,661,738,856,1486,1696,5192,5266,5732];
neuronIDType = [1, 3,  3,  1,  3,  3,  2,  2,  3,  1,  2,   1,   2,   2,   2,   1];

imageFormat = 'epsc';
outputFolder = '';
fh = 1;

%% 

figure(fh); clf; hold on;
colors = hsv(4);

nNeurons = length(neuronIDList);
for kk=1:nNeurons
    load(fullfile(statisticsActivationFolder, ...
        sprintf('neuron%d_act.mat',neuronIDList(kk))));

    % Finding where the data is stored
    stimDataPos = ismember(activationData.labels,'nSpikes+');
    stimDataNeg = ismember(activationData.labels,'nSpikes-');
    stimDataIrr = ismember(activationData.labels,'Irradiance-mW/mm^2');
    stimDataPW = ismember(activationData.labels,'Duration');
    
    % Keep only 4ms PW
    pwOfInterestPos = activationData.data(:,stimDataPW)==4;
    nSpikesPos = activationData.data(pwOfInterestPos,stimDataPos);
    nSpikesNeg = abs(activationData.data(pwOfInterestPos,stimDataNeg));
    irrValues = activationData.data(pwOfInterestPos,stimDataIrr);
    
    % Sort things
    [nSpikesPosSorted, ind] = sort(nSpikesPos);
    nSpikesNegSorted = nSpikesNeg(ind);

    % Plots
    neuronColor = colors(neuronIDType(kk), :);
    scatter(nSpikesPosSorted, nSpikesNegSorted, 'MarkerFaceColor',neuronColor,...
        'MarkerEdgeColor',neuronColor,'Marker','square');
    plot(nSpikesPosSorted, nSpikesNegSorted,'color',neuronColor);
end

line([0 5], [0 5], 'color', 'r')
axis equal
xlim([0 3])
ylim([0 3])
xlabel('Number of elicited action potentials');
ylabel('Number of action potentials suppressed');
