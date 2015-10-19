dataPath = '/Volumes/Analysis/2012-09-24-3/data007/data007'; 

elecRespPath = '/Volumes/Analysis/2012-09-24-3/443_elec_resp/'; 

autoFilePath = '/Volumes/Analysis/2012-09-24-3/data006-autosort/';

datarun = load_data(dataPath);

datarun = load_neurons(datarun);

datarun = load_sta(datarun, 'load_sta', 'all');

datarun = load_params(datarun);

datarun = load_ei(datarun, 'all');

cellIDs = showCellsNearElec(datarun,443,elecRespPath);

showCellsNearElec(datarun,443,autoFilePath);
