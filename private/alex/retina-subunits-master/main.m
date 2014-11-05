%% cone preprocessing

datarun = load_data(fullfile(server_path(), '2011-10-25-5/streamed/data001-0/data001-0'));
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun = load_cones(datarun,'15'); % load_cones(datarun, 'Analysis');

datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
conepreprocess_save(datarun, 'cone_data_ind', 'data','date','2011-10-25-5');


%% cone preprocessing

%datarun = load_data(fullfile(server_path(), '2012-09-24-1/data003/data003'));
%datarun = load_params(datarun,'verbose',1);
%datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
%datarun = set_polarities(datarun);
%datarun = load_neurons(datarun);

%datarun = load_cones(datarun,'15');

%datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
%conepreprocess_save(datarun, 'cone_data_ind', 'data','date','2012-09-24-1');
    
%% cone preprocessing

datarun = load_data(fullfile(server_path(), '2011-06-30-0/streamed/data003/data003'));
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun = load_cones(datarun,'10');

datarun.names.rrs_movie_path='/Volumes/Analysis/2011-06-30-0/data003/data003.movie';%BW-2-4-0.48-11111-300x300-60.35.xml';

datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
conepreprocess_save(datarun, 'cone_data_ind', 'data', 'date', '2011-06-30-0');

    
%%
dat = loadData('bertha','2011-06-30-0/data003');
batchAnal('bertha','2011-06-30-0/data003','off midget','stc',1,0,0.33);
batchAnal('bertha','2011-06-30-0/data003','off midget','subunit',1,0,0.33)
batchAnal('bertha','2011-06-30-0/data003','off midget','subunit',0,1,0.33)


%%
%dat = loadData('stanford','2011-10-25-5/data001-0');
%batchAnal('stanford','2011-10-25-5/data001-0','off midget','stc',1,0,0.33);
%batchAnal('stanford','2011-10-25-5/data001-0','off midget','subunit',1,0,0.33)
%batchAnal('stanford','2011-10-25-5/data001-0','off midget','subunit',0,1,0.33)

%%
%2012-09-24-1/data003

workstation='stanford';
% dataset='2013-08-19-2/data001';
dataset='2011-10-25-5/data001-0';

dat = loadData(workstation,dataset);
disp('now running stc')
batchAnal(workstation,dataset,'off midget','stc',1,0,0.33);
disp('now running subunits')
tic
batchAnal(workstation,dataset,'off midget','subunit-local',1,0,0.33)
toc
disp('now running figures')
batchAnal(workstation,dataset,'off midget','subunit-local',0,1,0.33)


