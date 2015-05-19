%% load stuff
datarun = load_data('/Volumes/Analysis/2010-09-24-1/data006/data006');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun1 = load_data('/Volumes/Analysis/2010-09-24-1/data021/data021');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);
datarun1 = load_neurons(datarun1);


datarun2 = load_data('/Volumes/Analysis/2010-09-24-1/data035/data035');
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2);
datarun2 = set_polarities(datarun2);
datarun2 = load_neurons(datarun2);

datarun3 = load_data('/Volumes/Analysis/2010-09-24-1/data002/data002');
datarun3 = load_params(datarun3,'verbose',1);
datarun3 = load_sta(datarun3);
datarun3 = set_polarities(datarun3);
datarun3 = load_neurons(datarun3);

datarun4 = load_data('/Volumes/Analysis/2010-09-24-1/data005/data005');
datarun4 = load_params(datarun4,'verbose',1);
datarun4 = load_sta(datarun4);
datarun4 = set_polarities(datarun4);
datarun4 = load_neurons(datarun4);

%% find stable cells

figure
plot_rf_summaries(datarun, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun1, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')
plot_rf_summaries(datarun2, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'g')
plot_rf_summaries(datarun3, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'k')

figure
plot_rf_summaries(datarun3, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun4, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')


figure
plot_rf_summaries(datarun, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun4, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')


figure
plot_rf_summaries(datarun, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun2, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')


figure
plot_rf_summaries(datarun1, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun2, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')



%% load stuff
datarun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data006-from-d05_36/data006-from-d05_36');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun1 = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data035-from-d05_36/data035-from-d05_36');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);
datarun1 = load_neurons(datarun1);

vormap = load('/Volumes/Archive/2010-09-24-1/Visual/cone maps data006/full_cone_map.txt');
figure
colormap gray
imagesc(vormap)


vorrun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data034-from-d05_36/data034-from-d05_36');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
[inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-gaussian-1-6-0.16-11111-2023x1.xml');

%% find stable runs

tmp_map = vormap;
for i=1:max(tmp_map(:))  
    tmp_map(tmp_map==i) = 0.3+(rand(1)-0.5)/8;
end

cell_type = 4;
bords = 10;
figure
set(gcf, 'position', [58 102 1182 996])
cnt=1;
tmp = zeros(320,320,3);
tmp(:,:,2) = tmp_map;

for i=109:117%find(ww./pixn<5)
    visionID = datarun.cell_types{cell_type}.cell_ids(i);
    datarunID = find(datarun.cell_ids == visionID);
    sta = double(datarun.stas.stas{datarunID}(:,:,1,4));
    sta1 = double(datarun1.stas.stas{datarunID}(:,:,1,4));
    
%     figure
%     colormap gray
%     imagesc(sta)
    
    if datarun.stas.polarities{datarunID} == 1
        [tt, ind] = max(sta(:));
        [tt1, ~] = max(sta1(:));
    else
        [tt, ind] = min(sta(:));
        [tt1, ~] = min(sta1(:));
    end
    [a,b] = ind2sub([320,320],ind);
    
    sta = sta/tt;
    sta1 = sta1/tt1;

    tmp(:,:,1) = sta;
    tmp(:,:,3) = sta1; 
    figure(1)
    subplot(3,3,cnt)
    imagesc(tmp)
    axis([b-bords b+bords a-bords a+bords])    
    title(['cell ID ', int2str(datarunID), ', vision ID ', int2str(visionID), ', i=', int2str(i)])
    
%     tt = unique(vormap(find(sta>0.4)));
%     tt(tt==0)=[];
%     clear sta_inds sta1_inds
%     for j = 1:length(tt)
%         a = find(vormap == tt(j));
%         tmp1 = sta(a);
%         tmp1 = tmp1(tmp1>0.4);
%         sta_inds(j)=sum(tmp1);
%         tmp1 = sta1(a);
%         tmp1 = tmp1(tmp1>0.4);
%         sta1_inds(j)=sum(tmp1);
%         [~, ic] = sort(sta_inds);
%     end
%     figure(2)
%     subplot(3,4,cnt)
%     plot(sta_inds(ic), sta1_inds(ic), '*')    
%     title(['cell ID ', int2str(cellInd), ', vision ID ', int2str(visionID)])    
    
    cnt = cnt+1;
end


good_cells = [23, 64, 87, 100, 106, 122, 140, 145, 146, 149, 150, 153, 156, ...
    159, 164, 169, 170, 172, 174, 175, 176, 179, 180, 181, 182, 184, 187, ...
    191, 208];

%% biased STA
datarunID = 64;

visionID = datarun.cell_ids(datarunID);

sta_params.length = 6;
sta_params.offset = 0;

my_sta=zeros(size(inputs_v,1),sta_params.length);

spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];
nspikes = numel(spikes);

spike_rate_all=zeros(size(inputs_v,2),1);
    
while ~isempty(spikes)
    [~, ia, ~] = unique(spikes);
    spike_rate_all(spikes(ia))=spike_rate_all(spikes(ia))+1;
    for j=1:sta_params.length
        my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
            sum(inputs_v(:,spikes(ia)-sta_params.length+j+sta_params.offset),2);
    end
    spikes(ia)=[];
end
my_sta=my_sta/nspikes;

% get voronoi map sta
tmp_map = vormap;
tt=0:max(tmp_map(:));
vorsta=zeros(320,320,sta_params.length);
coneX = zeros(max(tmp_map(:)),1);
coneY = coneX;
for i=1:2023
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
for i=1:2023
    text(coneY(i),coneX(i), int2str(i), 'color', 'r')
end


% plot voronoi sta frame by frame
figure
for i=1:6
    subplot(2,3,7-i)
    colormap gray
    imagesc(vorsta(:,:,i))
    title(['frame ',int2str(i)])
end
