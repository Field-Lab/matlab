function [elecs] = getBundleElecs_overtime(rawData, movieNo, threshold, t)
%%Takes EI loaded by generateEiFromStimPattern for paricular pattern and checks which electrodes are involves at a given amplitude

%Average over trials and take min of time
rawData = squeeze(mean(rawData,1)); %average over trials
rawData = squeeze(rawData(:, t, movieNo));

%Check which EIs are below threshold at any point in time
elecs = find(rawData < threshold);
