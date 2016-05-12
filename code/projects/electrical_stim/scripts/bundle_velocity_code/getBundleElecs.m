function [elecs] = getBundleElecs(rawData, movieNo, threshold)
%%Takes EI loaded by generateEiFromStimPattern for paricular pattern and checks which electrodes are involves at a given amplitude

%Average over trials and take min of time
rawData = squeeze(mean(rawData,1));
rawData = squeeze(min(rawData,[],2));
rawData = rawData(:, movieNo); %now a single vector with 512 min values at each elec

%Check which EIs are below threshold at any point in time
elecs = find(rawData < threshold);
