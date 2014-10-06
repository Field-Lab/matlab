function datarun = createDatarunVar()
pathname = uigetdir('/Volumes/Analysis/','select a directory with vision output files'); 
dataPath = [pathname pathname(find(pathname == filesep,1,'last'):end)];
datarun  = load_data(dataPath);
datarun  = load_neurons(datarun);
datarun  = load_sta(datarun, 'load_sta', 'all');
datarun  = load_params(datarun);
datarun  = load_ei(datarun, 'all');
end