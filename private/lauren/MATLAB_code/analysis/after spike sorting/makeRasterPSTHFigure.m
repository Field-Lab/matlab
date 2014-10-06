
clear all

elecRespPath = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002';
neuronNo = 902;
patternNo = 62;
movieNo = 26;

rasterOptions.markerSize = 5;
rasterOptions.markerColor = [0 0 0];

psthOptions.binSize = 0.1; % in ms
psthOptions.lineColor = [0 0 0];

% options.markerSize = 5; %for raster plot
% options.markerColor = [0 0 0];
% options.binSize = 0.1; %in ms
% options.lineColor = [0 0 0];



load([elecRespPath '/elecResp_n' num2str(neuronNo) '_p'...
    num2str(patternNo) '.mat'])

figure('position', [100 100 200 400])


axes('position', [0.18 0.55 0.8 0.4])
title(['neuron ' num2str(neuronNo) ', pattern ' num2str(patternNo)])
plot_raster(gca, elecResp, movieNo, rasterOptions)
xlabel('')


axes('position', [0.18 0.1 0.8 0.4])
plot_psth(gca, elecResp, movieNo, psthOptions)

