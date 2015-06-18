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

%% pre-calculate and save stuff out

% RGB gaussian voronoi

vorrun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data034-from-d05_36/data034-from-d05_36');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
[inputs_v, refresh_v, duration_v] = get_wn_movie_ath_rgb(vorrun, 'RGB-gaussian-1-6-0.16-11111-2023x1.xml');

all_stas = zeros(2023,3,6,length(vorrun.cell_ids));

sta_params.length = 6;
sta_params.offset = 0;

for datarunID = 1:length(vorrun.cell_ids)    
    
    my_sta=zeros(size(inputs_v,1),3,sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    nspikes = numel(spikes);
    
    while ~isempty(spikes)
        [~, ia, ~] = unique(spikes);
        for j=1:sta_params.length
            my_sta(:,:,sta_params.length-j+1) = my_sta(:, :,sta_params.length-j+1) +...
                sum(inputs_v(:,:,spikes(ia)-sta_params.length+j+sta_params.offset),3);
        end
        spikes(ia)=[];
    end
    my_sta=my_sta/nspikes;
    
    all_stas(:,:,:, datarunID) = my_sta;
end

save('/Users/alexth/Desktop/voronoi/2010-09-24-1/data034_stas.mat', 'all_stas')



% BW binary voronoi

vorrun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data036-from-d05_36/data036-from-d05_36');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
[inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-1-6-0.48-11111-2023x1.xml');

all_stas = zeros(2023,6,length(vorrun.cell_ids));

sta_params.length = 6;
sta_params.offset = 0;

for datarunID = 1:length(vorrun.cell_ids)    
    
    my_sta=zeros(size(inputs_v,1),sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    nspikes = numel(spikes);
    
    while ~isempty(spikes)
        [~, ia, ~] = unique(spikes);
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(:,spikes(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes(ia)=[];
    end
    my_sta=my_sta/nspikes;
    
    all_stas(:,:,datarunID) = my_sta;
end

save('/Users/alexth/Desktop/voronoi/2010-09-24-1/data036_stas.mat', 'all_stas')



%% load stuff
datarun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data006-from-d05_36/data006-from-d05_36');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

% [inputs, refresh, duration] = get_wn_movie_ath_rgb(datarun, 'BW-1-6-0.48-11111-320x320.xml');

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
% vorrun = load_sta(vorrun);
[inputs_v, refresh_v, duration_v] = get_wn_movie_ath_rgb(vorrun, 'RGB-gaussian-1-6-0.16-11111-2023x1.xml');


vorrun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data036-from-d05_36/data036-from-d05_36');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
vorrun = load_sta(vorrun);
[inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-1-6-0.48-11111-2023x1.xml');
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
% 
% figure
% colormap gray
% imagesc(vorsta(:,:,3))
% hold on
% for i=1:2023
%     text(coneY(i),coneX(i), int2str(i), 'color', 'r')
% end


% plot voronoi sta frame by frame
figure
for i=1:6
    subplot(2,3,7-i)
    colormap gray
    imagesc(vorsta(:,:,i))
    title(['frame ',int2str(i)])
end

figure
plot(my_sta')


%% Plot stuff
clear 
clc
% white noise, first run
datarun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data006-from-d05_36/data006-from-d05_36');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

mdf_file=['/Volumes/Analysis/stimuli/white-noise-xml/', 'BW-1-6-0.48-11111-320x320.xml'];
[~,~,~,duration,refresh] = get_movie_ath(mdf_file,datarun.triggers, 1,2);

% white nooise, second run
datarun1 = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data035-from-d05_36/data035-from-d05_36');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);
datarun1 = load_neurons(datarun1);

mdf_file=['/Volumes/Analysis/stimuli/white-noise-xml/', 'BW-1-6-0.48-22222-320x320.xml'];
[~,~,~,duration1,refresh1] = get_movie_ath(mdf_file,datarun1.triggers, 1,2);

% voronoi, RGB Gaussian, first run
vorrun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data034-from-d05_36/data034-from-d05_36');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
load('/Users/alexth/Desktop/voronoi/2010-09-24-1/data034_stas.mat', 'all_stas');
sta_34=all_stas;

mdf_file=['/Volumes/Analysis/stimuli/white-noise-xml/', 'RGB-gaussian-1-6-0.16-11111-2023x1.xml'];
[~,~,~,duration_v,refresh_v] = get_movie_ath(mdf_file,vorrun.triggers, 1,2);

% voronoi, BW binary, second run
vorrun1 = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data036-from-d05_36/data036-from-d05_36');
vorrun1 = load_params(vorrun1,'verbose',1);
vorrun1 = load_neurons(vorrun1);
load('/Users/alexth/Desktop/voronoi/2010-09-24-1/data036_stas.mat', 'all_stas');
sta_36=all_stas;

mdf_file=['/Volumes/Analysis/stimuli/white-noise-xml/', 'BW-1-6-0.48-11111-2023x1.xml'];
[~,~,~,duration_v1,refresh_v1] = get_movie_ath(mdf_file,vorrun1.triggers, 1,2);

% voronoi general map
vormap = load('/Volumes/Archive/2010-09-24-1/Visual/cone maps data006/full_cone_map.txt');
tmp_map = vormap;
for i=1:max(tmp_map(:))  
    tmp_map(tmp_map==i) = 0.3+(rand(1)-0.5)/8;
end

clear mdf_file all_stas

% use temporal summation for highest quality

filepath = '/Users/alexth/Desktop/voronoi/2010-09-24-1/d05-36-norefit/';
borders(1) = 20;
borders(2) = 20;
borders(3) = 10;
borders(4) = 10;
borders(5:15) = 30;
for datarunID = 1:length(datarun.cell_ids)
    visionID = datarun.cell_ids(datarunID);
    [folder, my_type] = find_cell_type(datarun, visionID);
    
    if ~exist([filepath,folder],'dir')
        mkdir([filepath,folder]);
    end
    
    
    fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[1 114 1618 984]);
    
    % first white noise
    subplot(2,3,1)
    sta = double(squeeze(datarun.stas.stas{datarunID})); 
    tmp = reshape(sta(:,:,1:4), 320*320,4);
    if datarun.stas.polarities{datarunID}<0
        tmp = -tmp;
        sign_sta = -1;
    else
        sign_sta = 1;
    end    
    thr = robust_std(tmp(:,4))*2;
    thr1 = -robust_std(tmp(:,2))*2;      
    tc = tmp(tmp(:,2)<thr1 & tmp(:,4)>thr,:);
    if numel(tc)>4
        tc = mean(tc);
    end
    if isempty(tc)
        tmp = tmp(:,4);
        tc_leg = 'slice';
    else
        tmp = tmp.*repmat(tc, size(tmp,1),1);
        tmp = sum(tmp,2);
        tc_leg = 'full';
    end
    tmp = sign_sta*reshape(tmp, 320,320);
    
    colormap gray
    imagesc(tmp)
    comp_sta = tmp;
    title(['first white noise, data006, tc ', tc_leg])
    
    % second white noise
    subplot(2,3,4)
    sta = double(squeeze(datarun1.stas.stas{datarunID})); 
    tmp = reshape(sta(:,:,1:4), 320*320,4);
    if datarun1.stas.polarities{datarunID}<0
        tmp = -tmp;
        sign_sta = -1;
    else
        sign_sta = 1;
    end    
    thr = robust_std(tmp(:,4))*2;
    thr1 = -robust_std(tmp(:,2))*2;
    tc = tmp(tmp(:,2)<thr1 & tmp(:,4)>thr,:);
    if numel(tc)>4
        tc = mean(tc);
    end
    if isempty(tc)
        tmp = tmp(:,4);
        tc_leg = 'slice';
    else
        tmp = tmp.*repmat(tc, size(tmp,1),1);
        tmp = sum(tmp,2);
        tc_leg = 'full';
    end
    tmp = sign_sta*reshape(tmp, 320,320);
    
    colormap gray
    imagesc(tmp)
    comp_sta1 = tmp;
    title(['second white noise, data035, tc ', tc_leg])
    
    
    % first voronoi (Gaussian RGB)
    subplot(2,3,2)
    sta = sta_34(:,:,:,datarunID);     
    
    if datarun.stas.polarities{datarunID}<0
        sta = -sta;
        sign_sta = -1;
    else
        sign_sta = 1;
    end    
    % time courses separately for color channels
    color_sta = zeros(size(sta_34,1),3);
    get_single = 1;
    if ~isempty(datarun.stas.polarities{datarunID})
        for j = 1:3
            tmp = squeeze(sta(:,j,2:5));
            thr1 = -robust_std(tmp(:,3))*2;
            thr = robust_std(tmp(:,1))*2;
            tc = tmp(tmp(:,3)<thr1 & tmp(:,1)>thr,:);
            if numel(tc)>4
                tc = mean(tc);
            end
            if ~isempty(tc)
                get_single = 0;
                tc_leg = 'full';
                tmp = tmp.*repmat(tc, size(tmp,1),1);
                tmp = sum(tmp,2);
                color_sta(:,j) = tmp;
            else
                get_single = 1;
            end
        end
    end
    if get_single
        tc_leg = 'slice';
        for j = 1:3
            tmp = squeeze(sta(:,j,2));
            color_sta(:,j) = tmp;
        end
    end
    

    color_sta = color_sta*sign_sta;    
    abs_max =  0.5/max(abs(color_sta(:)));    
    % get voronoi map sta
    tt=0:max(vormap(:));
    vorsta=zeros(320,320,3);
    cnt = 1;
    for i=1:2023
        [a, b] = find(vormap==tt(i+1));
        if ~isempty(a)
            for j = 1:length(a)
                vorsta(a(j),b(j),:) = color_sta(cnt,:);
            end
        end
        cnt=cnt+1;
    end   
    imagesc(vorsta*abs_max+0.5)
    title(['first voronoi, data034 (RGB Gaussian), tc ', tc_leg])
    
    
    % second voronoi (Binary BW)
    subplot(2,3,5)
    sta = sta_36(:,:,datarunID);    
    if datarun.stas.polarities{datarunID}<0
        sta = -sta;
        sign_sta = -1;
    else
        sign_sta = 1;
    end
    get_single = 1;
    if ~isempty(datarun.stas.polarities{datarunID})
        sta = sta(:,2:5);
        thr = robust_std(sta(:,1))*2;
        thr1 = -robust_std(sta(:,3))*2;
        tc = sta(sta(:,3)<thr1 & sta(:,1)>thr,:); 
        if numel(tc)>4
            tc = mean(tc);
        end
        if ~isempty(tc)
            sta = sta.*repmat(tc, size(sta,1),1);
            sta = sum(sta,2);
            sta = sign_sta*sta;
            tc_leg = 'full';
            get_single = 0;
        else
            get_single = 1;
        end
    end
    if get_single
        sta = sta(:,2);
        tc_leg = 'slice';
    end
    % get voronoi map sta
    tt=0:max(vormap(:));
    vorsta=zeros(320,320);
    cnt = 1;
    for i=1:2023
        [a, b] = find(vormap==tt(i+1));
        if ~isempty(a)
            for j = 1:length(a)
                vorsta(a(j),b(j)) = sta(cnt);
            end
        end
        cnt=cnt+1;
    end
    colormap gray
    imagesc(vorsta)
    title(['second voronoi, data036 (BW binary), tc ', tc_leg])  
    
    
    % combined plot for quality access    
    subplot(2,3,3)
    
    % find zone of interest    
    tmp = comp_sta+comp_sta1;
    tt = max(abs(tmp(:)));
    [a,b] = find(abs(tmp)==tt);
    a = max(1,a-borders(my_type)):min(320,a+borders(my_type));
    b = max(1,b-borders(my_type)):min(320,b+borders(my_type));
    
    % combine channels    
    combo_sta = zeros(length(a),length(b),3);
    combo_sta(:,:,1) = sign_sta*comp_sta(a,b)/max(abs(comp_sta(:)));
    combo_sta(:,:,3) = sign_sta*comp_sta1(a,b)/max(abs(comp_sta1(:)));
    combo_sta(:,:,2) = tmp_map(a,b);
    
    imagesc(combo_sta);
    
    title('white noise combo over voronoi map')
    
    
    % SNR in first and second WN
    subplot(2,3,6)
    tt = sort(comp_sta(:));
    plot(tt([1:200 end-200:end]), 'linewidth',2) 
    hold on
    tt = sort(comp_sta1(:));
    plot(tt([1:200 end-200:end]), 'linewidth',2)
    h = legend('first WN', 'second WN');
    set(h,'location', 'best')
    line([0 400], [0, 0], 'color','k')
    axis tight
    
    title([folder, ', cell ', int2str(visionID), ' (datarun ID ',int2str(datarunID), ')'])
    
    

    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(datarunID)]));
    close(fig)

end
