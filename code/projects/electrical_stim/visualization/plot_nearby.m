% dataPath = '/Volumes/Analysis/2012-09-24-3/data007/data007'; 
% elecRespPath = '/Volumes/Analysis/2012-09-24-3/443_elec_resp/'; 
% autoFilePath = '/Volumes/Analysis/2012-09-24-3/data006-autosort/';
dataPath = '/Volumes/Analysis/2015-11-09-3/data000/data000'; 
autoFilePath = '/Volumes/Analysis/2015-11-09-3/data001-data002-autosort/';
autoFilePath2 = '/Volumes/Analysis/2015-11-09-3/data003-autosort/';
datarun = load_data(dataPath);
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');

cellIDs = showCellsNearElec(datarun,443,'datapath',elecRespPath);
showCellsNearElec(datarun,443,'datapath',autoFilePath,'autoFile',true);
showCellsNearElec(datarun,443,'datapath',autoFilePath,'autoFile',true);

showCellsNearElec(datarun,320,'datapath',autoFilePath,'autoFile',true,'threshold',10);
showCellsNearElec(datarun,26,'datapath',autoFilePath,'autoFile',true,'threshold',10);
showCellsNearElec(datarun,83,'datapath',autoFilePath,'autoFile',true,'threshold',10);
showCellsNearElec(datarun,200,'datapath',autoFilePath,'autoFile',true,'threshold',10);
showCellsNearElec(datarun,443,'datapath',autoFilePath,'autoFile',true,'threshold',10);


dataPath = '/Volumes/Analysis/2015-10-06-3/data000/data000'; 
autoFilePath = '/Volumes/Analysis/2015-10-06-3/data001-data002-autosort/';
datarun = load_data(dataPath);
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');
