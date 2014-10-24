%% Normal white noise
datarun_wn = load_data('2012-08-09-1/streamed/data001/data001');
datarun_wn = load_params(datarun_wn);
datarun_wn = load_neurons(datarun_wn);
datarun_wn = load_cones(datarun_wn);
%datarun_wn = load_java_movie(datarun_wn);
datarun_wn = load_ei(datarun_wn,'all');

datarun_wn = conepreprocess_wnm(datarun);


%% White noise repeats
datarun = load_data('2012-08-09-1/data002');
datarun = load_params(datarun);
datarun = load_neurons(datarun);

moviexmlpath = '/share/users-lcv/freeman/retinalData/stimuli/BW-3-5-0.48-11111-200x200-60.35.xml';
triggers_per_trial = 6;
conepath = '2012-08-09-1_data001';

datarun = conepreprocess_wnrepeats(datarun, moviexmlpath, triggers_per_trial, 'conepath', conepath);
datarun = load_ei(datarun);