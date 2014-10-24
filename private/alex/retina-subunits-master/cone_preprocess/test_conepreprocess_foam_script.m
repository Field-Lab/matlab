%%
clear all
%% load the cone finding run
basePath = '/share/users-lcv/freeman/retinalData/';
datarun_wn = load_data(strcat(basePath,'rawdata/2012-08-21-0/streamed/data001/data001'));
datarun_wn = load_params(datarun_wn);
datarun_wn = load_cones(datarun_wn);

%% load the validation run (e.g. foam)
datarun = load_data(strcat(basePath,'rawdata/2012-08-21-0/data003/data003'));
datarun = load_params(datarun);
datarun.cones = datarun_wn.cones;
datarun.names.cones = datarun_wn.names.cones;
%%
datarun = conepreprocess_foam(datarun);
%%
datarun = conepreprocess_wnm(datarun_wn);