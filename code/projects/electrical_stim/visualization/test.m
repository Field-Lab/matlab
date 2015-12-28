dataPath = 'E:\2015-11-09-10\data001';

%Make datarun variable
datarun  = load_data(dataPath); %this isn't doing its job on Windows for some reason, so doing it manually - Sasi
dataPathSplit = strsplit(dataPath, filesep);
datarun.names.rrs_neurons_path = [dataPath filesep dataPathSplit{end} '.neurons'];
datarun.names.rrs_ei_path = [dataPath filesep dataPathSplit{end} '.ei'];
datarun  = load_neurons(datarun);
datarun  = load_ei(datarun, 'all');
%datarun.names.rrs_sta_path = 'E:\2015-11-09-10\data001\data001.sta';
%datarun.names.rrs_params_path = 'E:\2015-11-09-10\data001\data001.params';
%datarun  = load_sta(datarun, 'load_sta', 'all');
%datarun  = load_params(datarun);

%loop over thresholds
thr = -50; 
[cellIdsToCheck, cellIndices, electrodes] = getLargeAmpSpikes(datarun, thr); 

%generate AllEIs
numNeurons         = size(datarun.ei.eis,1); 
for n  = 1:numNeurons
    allEIs(n,:,:) = datarun.ei.eis{n}; 
end
