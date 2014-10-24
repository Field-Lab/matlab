%% Normal white noise
moviexmlpath_wn = '/share/users-lcv/freeman/retinalData/stimuli/BW-3-5-0.48-11111-200x200-60.35.xml';

datarun_wn = load_data('2012-09-18-3/streamed/data003/data003');
datarun_wn = load_params(datarun_wn);
datarun_wn = load_neurons(datarun_wn);
datarun_wn = load_cones(datarun_wn);
datarun_wn.names.movie_xml_path = moviexmlpath_wn;
datarun_wn = load_java_movie(datarun_wn);
datarun_wn = load_ei(datarun_wn,'all');

datarun_wn = conepreprocess_wnm(datarun);


%% White noise repeats
datarun = load_data('2012-09-18-3/streamed/data004/data004');
datarun = load_params(datarun);
datarun = load_neurons(datarun);

moviexmlpath = '/share/users-lcv/freeman/retinalData/stimuli/BW-3-5-0.48-11111-200x200-60.35.xml';
triggers_per_trial = 6;
conepath = '2012-08-09-1_data001';

datarun = conepreprocess_wnrepeats(datarun, moviexmlpath, triggers_per_trial, 'conepath', conepath);
datarun = load_ei(datarun,'all');

cell_list_map = map_ei(datarun_wn,datarun);

datarun = rmfield(datarun,'ei');
datarun.cell_list_map = cell_list_map;

save('/share/users-lcv/freeman/retinalData/data/2012-08-09-1/data002/conepreprocess.mat','datarun')

