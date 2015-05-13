%% load stuff
datarun = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data008-from-d08_11/data008-from-d08_11');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);

datarun1 = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data011-from-d08_11/data011-from-d08_11');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);

vormap = load('/Volumes/Data/2011-12-13-2/Visual/2011-12-13-2_f04_vorcones/map-0000.txt');
figure
imagesc(vormap)


vorrun = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data009-from-d08_11/data009-from-d08_11');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
[inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-937x1-60.35.xml');

%% get STA


my_cell = 27;

sta_params.length = 15;
sta_params.offset = 0;

my_sta=zeros(size(inputs,1),sta_params.length);

cellInd = find(datarun.cell_ids == datarun.cell_types{4}.cell_ids(my_cell));
spikes=ceil((vorrun.spikes{cellInd}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];
nspikes = numel(spikes);

spike_rate_all=zeros(size(inputs,2),1);
    
while ~isempty(spikes)
    [~, ia, ~] = unique(spikes);
    spike_rate_all(spikes(ia))=spike_rate_all(spikes(ia))+1;
    for j=1:sta_params.length
        my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
            sum(inputs(:,spikes(ia)-sta_params.length+j+sta_params.offset),2);
    end
    spikes(ia)=[];
end
my_sta=my_sta/nspikes;



% get voronoi map sta
tmp_map = vormap;
tt=0:max(tmp_map(:));
vorsta=zeros(600,600,sta_params.length);
coneX = zeros(max(tmp_map(:)),1);
coneY = coneX;
for i=1:937
    [a, b] = find(tmp_map==tt(i+1));
    % find center of this cone
    coneX(i) = mean(a);
    coneY(i) = mean(b);
    if ~isempty(a)
        for j = 1:length(a)
            vorsta(a(j),b(j),:) = my_sta(i,:);
        end
    end
end

figure
colormap gray
imagesc(vorsta(:,:,3))
hold on
for i=1:937
    text(coneY(i),coneX(i), int2str(i), 'color', 'r')
end



% plot voronoi sta frame by frame
figure
for i=1:9
    subplot(3,3,10-i)
    colormap gray
    imagesc(vorsta(:,:,i))
    title(['frame ',int2str(i)])
end



%plot sta to find the threshold
figure
plot(my_sta')
% cones with real center input
sel_cones.real = find(my_sta(:,3)<-0.041);
% cones with surround input
sel_cones.sur = find(my_sta(:,3)>0.017);
% cones with no iput (furthest away)
realX = mean(coneX(sel_cones.real));
realY = mean(coneY(sel_cones.real));
tmp = pdist2([realX, realY], [coneX coneY]);
[~, ic] = sort(tmp);
ic(isnan(tmp(ic))) = [];
sel_cones.far = ic(end-25:end)';

cones = [sel_cones.far; sel_cones.real; sel_cones.sur];
vorsta_select_cones=zeros(600,600);
for i=cones'
    [a, b] = find(tmp_map==tt(i+1));
    if ~isempty(a)
        for j = 1:length(a)
            vorsta_select_cones(a(j),b(j)) = my_sta(i,3);
        end
    end
end

figure
subplot(1,2,1)
plot(my_sta(cones,:)')
title([int2str(length(cones)), ' cones'])

subplot(1,2,2)
colormap gray
imagesc(vorsta_select_cones)
hold on
for i=cones'
    text(coneY(i),coneX(i), int2str(i), 'color', 'r')
end


tmprun = vorrun;

datarun1 = get_sta_summaries(datarun1, {1,2,3,4,5}, 'keep_stas', false, 'marks_params', struct('robust_std_method', 6));

my_filtered_inputs=zeros(size(inputs,1), length(spike_rate_all)-sta_params.length+1);  
for i=1:length(cones)
    my_filtered_inputs(i,:)=conv(inputs(cones(i),:),my_sta(i,:),'valid');
end

tmprun.stas.time_courses = cell(size(tmprun.spikes,1),1);
tmprun.stas.polarities = cell(size(tmprun.spikes,1),1);
tmprun.stas.rf_coms = cell(size(tmprun.spikes,1),1);
for i=1:size(tmprun.spikes,1)
    tmprun.stas.time_courses{i} = my_sta(i,end:-1:1);
    tmprun.stas.polarities{i} = datarun1.stas.polarities{i};
    tmprun.stas.rf_coms{i} = datarun1.stas.rf_coms{i}*2;
end
tmprun.spike_rate = zeros(size(tmprun.spikes,1),length(spike_rate_all)-sta_params.length+1);
tmprun.spike_rate(my_cell,:) = spike_rate_all(sta_params.length:end);
tmprun.cone_inputs = my_filtered_inputs';
datarun = tmprun;
save('/Volumes/Analysis/2011-12-13-2/subunits/data009/conepreprocess.mat','datarun')

%% run model


dataset='2011-12-13-2/data009';
rgcs=168;
cellTypeName='OFF midget';

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



