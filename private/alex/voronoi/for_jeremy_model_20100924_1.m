
%% create normal cone_preprocces on single cone white noise
datarun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data006-from-d05_36/data006-from-d05_36');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_cones_ath(datarun,'norefit');

datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'd06');
conepreprocess_save(datarun, 'cone_data_ind', 'bayes');

%% create cone_preprocess based on voronoi

% load normal WN run and get java sta fits
load('/Volumes/Analysis/2010-09-24-1/subunits/data006-from-d05_36/conepreprocess.mat')
datarun = get_sta_summaries(datarun, {1,2,3,4,5}, 'keep_stas', false, 'marks_params', struct('robust_std_method', 6));

my_cells = [827 1157 1907 2341 2851 3346 3961 4081 4306 4460 4637 4681 4862 ...
    4876 4936 4998 5011 5176 5191 5251 5476 5701 5851 5911 6076 6151 6241 6332 7036 ...
    4562 4611 5177 5206 5761 6511 5992 3202];

vorrun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data036-from-d05_36/data036-from-d05_36');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
wn_movie_name = 'BW-1-6-0.48-11111-2023x1.xml';
[inputs_v, refresh_v, duration_v] = get_wn_movie_ath(datarun, wn_movie_name);

vormap = load('/Volumes/Archive/2010-09-24-1/Visual/cone maps data006/full_cone_map.txt');

% load pre-calculated sta
load('/Users/alexth/Desktop/voronoi/2010-09-24-1/data036_stas.mat', 'all_stas');
sta_36=all_stas;

% temporary datarun with modified fields
modelrun = vorrun;

modelrun.stas.time_courses = cell(size(modelrun.spikes,1),1);
modelrun.stas.polarities = cell(size(modelrun.spikes,1),1);
modelrun.stas.rf_coms = cell(size(modelrun.spikes,1),1);
modelrun.stas.rfs = cell(size(modelrun.spikes,1),1);
modelrun.cones.weights = zeros(size(inputs_v,1),size(modelrun.spikes,1));
for i=1:size(modelrun.spikes,1)    
    modelrun.stas.polarities{i} = datarun.stas.polarities{i};
    modelrun.stas.rf_coms{i} = datarun.stas.rf_coms{i}*2;    
end
modelrun.cone_inputs = inputs_v';
modelrun.spike_rate = zeros(size(modelrun.spikes,1),size(inputs_v,2));


for i=1:length(my_cells)
    datarunID = find(datarun.cell_ids==my_cells(i));
    
    % calculate spike rate in frames
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spikes(spikes>size(inputs_v,2)) = []; % remove extra spikes (how the hell did they get there??)
    spikes(spikes<1) = [];
    spike_rate=zeros(size(inputs_v,2),1);    
    while ~isempty(spikes)
        [~, ia, ~] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    modelrun.spike_rate(datarunID,:) = spike_rate;
    
    % calculate time course from voronoi sta
    sta = all_stas(:,:,datarunID);   
    sign_sta = datarun.stas.polarities{datarunID};
    sta = sta*sign_sta;
    thr = robust_std(sta(:,2))*2;
    thr1 = -robust_std(sta(:,4))*2;
    tc = sta(sta(:,2)>thr & sta(:,4)<thr1,:);
    if numel(tc)>4
        tc = mean(tc);
    end
    modelrun.stas.time_courses{datarunID} = tc(end:-1:1)'* sign_sta; 
    
    % calculate cone weights (inner product of sta with time course); sta is
    % ON-like, so weights are positive for central cones
    sta = sta.*repmat(tc, size(sta,1),1);
    sta = sum(sta,2);
    modelrun.cones.weights(:,datarunID) = sta;
    
    % calculate rf (full voronoi map)
    tt=0:max(vormap(:));
    fullsta=zeros(320,320);
    cnt = 1;
    for k=1:2023
        [a, b] = find(vormap==tt(k+1));
        if ~isempty(a)
            for j = 1:length(a)
                fullsta(a(j),b(j)) = sta(cnt);
            end
        end
        cnt=cnt+1;
    end
    modelrun.stas.rfs{datarunID} = fullsta;    
    
end

% find cone centers
coneX = zeros(max(vormap(:)),1);
coneY = coneX;
for i=1:2023
    [a, b] = find(vormap==tt(i+1));
    coneX(i) = mean(a);
    coneY(i) = mean(b);
end
modelrun.cones.centers = [coneY coneX];
modelrun.cones.types(1:2023) = 'U';
modelrun.cones.types = modelrun.cones.types';

% save 
datarun = modelrun;
savepath = '/Volumes/Analysis/2010-09-24-1/subunits/data036-from-d05_36/';
mkdir(savepath);
save([savepath,'conepreprocess.mat'],'datarun')


%% shorten run

load('/Volumes/Analysis/2010-09-24-1/subunits/data036-from-d05_36/conepreprocess.mat')
time_border = 1500;
datarun.cone_inputs = datarun.cone_inputs(1:time_border,:);
datarun.spike_rate = datarun.spike_rate(:,1:time_border);
save(['/Volumes/Analysis/2010-09-24-1/subunits/data036-from-d05_36_to_',int2str(time_border),'/conepreprocess.mat'], 'datarun')

%% run model


dataset='2010-09-24-1/data036-from-d05_36';

rgcs = [827 1157 1907 2341 2851 3346 3961 4081 4306 4460 4637 4681 4862 ...
    4876 4936 4998 5011 5176 5191 5251 5476 5701 5851 5911 6076 6151 6241 6332 7036];

cellTypeName = 'OFF midget';
perc = 0.5;
workstation='bertha';

dat = loadData(workstation,dataset);

disp('now running stc')
batchAnal(workstation,dataset,cellTypeName,'stc',1,0,perc, rgcs);

disp('now running subunits')
tic
batchAnal(workstation,dataset,cellTypeName,'subunit-local',1,0,perc, rgcs)
toc

disp('now running figures')
batchAnal(workstation,dataset,cellTypeName,'subunit-local',0,1,perc, rgcs)



