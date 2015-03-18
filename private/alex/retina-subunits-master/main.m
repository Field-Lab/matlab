%% cone preprocessing

datarun = load_data(fullfile(server_path(), '2011-10-25-5/streamed/data001-0/data001-0'));
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun = load_cones_ath(datarun,'15'); % load_cones(datarun, 'Analysis');

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

datarun = load_cones_ath(datarun,'2011-06-30-0');

datarun.names.rrs_movie_path='/Volumes/Analysis/2011-06-30-0/data003/data003.movie';%BW-2-4-0.48-11111-300x300-60.35.xml';

datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
conepreprocess_save(datarun, 'cone_data_ind', 'data', 'date', '2011-06-30-0');

    
%% cone preprocessing

datarun = load_data(fullfile(server_path(), '2012-09-13-2/data009/data009'));
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun = load_cones_ath(datarun,'localmax');

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

dataset='2012-04-13-1/data002'; % cellType 4, 'OFF midget'
rgcs=2389;
cellTypeName='OFF midget';

dataset='2012-09-06-0/data004'; %cellType 10,  'OFF sparsepicks'
rgcs=[5103 5746 6031 6317 7022];
cellTypeName='OFF sparsepicks';

dataset='2012-08-21-2/data001'; %cellType 4,  'OFF midget'
rgcs=[391 2266 2828 5447 7682];
cellTypeName='OFF midget';

dataset='2012-04-13-1/data002'; % cellType 4, 'OFF midget'
rgcs=7158;
cellTypeName='OFF midget';

dataset='2012-09-13-2/data009'; % cellType 5, 'SBC'
rgcs=[1321,1623,1982,2431,3138,7232,3035,5224,6271,5733,5137];
cellTypeName='SBC';

dataset='2012-09-21-2/data009'; % cellType 5, 'SBC'
rgcs=[2030 2315 2613 2630 3661];
cellTypeName='SBC';

dataset='2010-03-05-2/data001'; % cellType 5, 'SBC'
rgcs=[226,1576,2582,456,2326,3587,4352,4503,5270,4053];
cellTypeName='SBC';


dataset='2013-08-19-4/data001'; % cellType 5, 'SBC'
rgcs=[49,1488,2191,2671,3032,6842,5431,4321];
rgcs=[3032]
cellTypeName='SBC';


dataset = '2011-06-30-0/data003';
rgcs=[4218 3721 2237];
rgcs=2237;
cellTypeName='ON midget';

dataset = '2010-03-05-2/data013';
rgcs=[481,560,841,677];
rgcs=481;
cellTypeName='ON midget';

dataset = '2008-08-27-5/data003';
rgcs=[1951, 2986];
cellTypeName='ON midget';

dataset = '2011-12-13-2/data000-0';
rgcs=[2658, 2659, 3856, 3061, 2251, 2266];
rgcs=[2658 2659 3856];
rgcs=3856;
cellTypeName='ON midget';


dataset='2011-10-25-5/data001-0';
rgcs=[6034];
cellTypeName='OFF parasol';

dataset='2010-09-24-1/data006';
rgcs=[1351 1037 1561];
cellTypeName='OFF parasol';

% workstation='stanford';
workstation='bertha';
dat = loadData(workstation,dataset);

disp('now running stc')
batchAnal(workstation,dataset,cellTypeName,'stc',1,0,0.5, rgcs);

disp('now running subunits')
tic
batchAnal(workstation,dataset,cellTypeName,'subunit-local',1,0,0.5, rgcs)
toc

disp('now running figures')
batchAnal(workstation,dataset,cellTypeName,'subunit-local',0,1,0.5, rgcs)


dataset='2011-06-30-0/data003';
cellTypeName='ON parasol';


workstation='bertha';
dat = loadData(workstation,dataset);
disp('now running stc')
batchAnal(workstation,dataset,cellTypeName,'stc',1,0,0.33);
disp('now running subunits')
tic
batchAnal(workstation,dataset,cellTypeName,'subunit-local',1,0,0.33)
toc
disp('now running figures')
batchAnal(workstation,dataset,cellTypeName,'subunit-local',0,1,0.33)