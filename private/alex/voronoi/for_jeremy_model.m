%% load stuff
datarun1 = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data008-from-d08_11/data008-from-d08_11');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);
datarun1 = load_neurons(datarun1);

% wn_movie_name = 'BW-2-6-0.48-11111-300x300-60.35.xml';
% [inputs, refresh, duration] = get_wn_movie_ath(datarun, wn_movie_name);

vormap = load('/Volumes/Data/2011-12-13-2/Visual/2011-12-13-2_f04_vorcones/map-0000.txt');

vorrun = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data009-from-d08_11/data009-from-d08_11');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
[inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-937x1-60.35.xml');

%% get STA

sta_params.length = 12;
sta_params.offset = 0;
time_border = 4500;

for datarunID = 168
    visionID = datarun1.cell_ids(datarunID);
    
    % full voronoi    
    my_sta=zeros(size(inputs_v,1),sta_params.length);    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(time_border,1);
    spikes_tmp = spikes(spikes>sta_params.length-sta_params.offset & spikes<time_border);
    full_voronoi_sta = zeros(600,600, sta_params.length);
    nspikes = length(spikes_tmp);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(:,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    my_sta=my_sta/nspikes;
    
    % get voronoi map sta
    tmp_map = vormap;
    tt=0:max(tmp_map(:));
    coneX = zeros(max(tmp_map(:)),1);
    coneY = coneX;
    for i=1:937
        [a, b] = find(tmp_map==tt(i+1));
        if ~isempty(a)
            coneX(i) = mean(a);
            coneY(i) = mean(b);
            for j = 1:length(a)
                full_voronoi_sta(a(j),b(j),:) = my_sta(i,:);
            end
        end
    end
end



% plot voronoi sta frame by frame
figure
for i=1:12
    subplot(3,4,13-i)
    colormap gray
    imagesc(full_voronoi_sta(:,:,i))
    title(['frame ',int2str(i)])
end

plot(my_sta(:,3))
inds = find(my_sta(:,3)<-0.05);

modelrun = vorrun;

datarun1 = get_sta_summaries(datarun1, {1,2,3,4,5}, 'keep_stas', false, 'marks_params', struct('robust_std_method', 6));

load('/Volumes/Analysis/2011-12-13-2/subunits/data008-from-d08_11/conepreprocess.mat')

modelrun.stas.time_courses = cell(size(modelrun.spikes,1),1);
modelrun.stas.time_courses{datarunID} = mean(my_sta(inds,end:-1:1))';
modelrun.stas.polarities = cell(size(modelrun.spikes,1),1);
modelrun.stas.rf_coms = cell(size(modelrun.spikes,1),1);
modelrun.stas.rfs = cell(size(modelrun.spikes,1),1);
modelrun.cones.weights = zeros(size(inputs_v,1),size(modelrun.spikes,1));
for i=1:size(modelrun.spikes,1)    
    modelrun.stas.polarities{i} = datarun1.stas.polarities{i};
    modelrun.stas.rf_coms{i} = datarun1.stas.rf_coms{i}*2;    
end
modelrun.stas.rfs{datarunID} = full_voronoi_sta(:,:,3)/100;
modelrun.spike_rate = zeros(size(modelrun.spikes,1),length(spike_rate));
modelrun.spike_rate(datarunID,:) = spike_rate;
modelrun.cone_inputs = inputs_v(:,1:time_border)';

plot(modelrun.stas.time_courses{datarunID})

% modelrun.spike_rate = zeros(size(modelrun.spikes,1),length(spike_rate)-sta_params.length+1);
% modelrun.spike_rate(datarunID,:) = spike_rate(sta_params.length:end);
% modelrun.cone_inputs = my_filtered_inputs';

modelrun.cones.centers = [coneY coneX];
modelrun.cones.types(1:937) = 'U';
modelrun.cones.types = modelrun.cones.types';
modelrun.cones.weights(:,datarunID) = -my_sta(:,3);

datarun = modelrun;
save('/Volumes/Analysis/2011-12-13-2/subunits/data009_4500frames/conepreprocess.mat','datarun')

%% short regular run
load('/Volumes/Analysis/2011-12-13-2/subunits/data008/conepreprocess.mat')
time_border = 1500;
datarun.cone_inputs = datarun.cone_inputs(1:time_border,:);
datarun.spike_rate = datarun.spike_rate(:,1:time_border);
save('/Volumes/Analysis/2011-12-13-2/subunits/data008_1500frames/conepreprocess.mat', 'datarun')
%% run model


dataset='2011-12-13-2/data008-from-d08_11_2';
dataset='2011-12-13-2/data008_2';
dataset='2011-12-13-2/data008_5min';
dataset='2011-12-13-2/data008_1500frames';

dataset = '2011-12-13-2/data009_short';
dataset = '2011-12-13-2/data009_4500frames';
rgcs=[3736]%datarun.cell_types{1, 4}.cell_ids  %[1321 2971 3136 3138 3736];
cellTypeName='OFF midget';
perc = 0.5;
% workstation='stanford';
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



