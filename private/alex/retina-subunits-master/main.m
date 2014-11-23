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

% workstation='stanford';
workstation='bertha';
% dataset='2013-08-19-2/data001';
% dataset='2011-10-25-5/data001-0';

dataset='2012-09-21-2/data009';
rgcs=[5086  6346 6406 7696]; % cellType 8, 'OFF sparse picks'
cellTypeName='OFF sparse picks';

dataset='2011-12-13-2/data008-0';
rgcs=[1321 1351 2251 3586 4576 5162]; %cellType 4, 'OFF midget'
cellTypeName='OFF midget';

dataset='2012-04-13-1/data006'; % cellType 4, 'OFF midget'
rgcs=6136;
cellTypeName='OFF midget';

dataset='2012-09-06-0/data004'; %cellType 10,  'OFF sparsepicks'
rgcs=[5103 5746 6031 6317 7022];
cellTypeName='OFF sparsepicks';

dataset='2012-08-21-2/data001'; %cellType 4,  'OFF midget'
rgcs=[391 2266 2828 5447 7682];
cellTypeName='OFF midget';


dat = loadData(workstation,dataset);
disp('now running stc')
batchAnal(workstation,dataset,cellTypeName,'stc',1,0,0.33, rgcs);
disp('now running subunits')
tic
batchAnal(workstation,dataset,cellTypeName,'subunit-local',1,0,0.33, rgcs)
toc
disp('now running figures')
batchAnal(workstation,dataset,cellTypeName,'subunit-local',0,1,0.33, rgcs)

