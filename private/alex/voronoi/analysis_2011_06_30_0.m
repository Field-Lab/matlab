%% load stuff
datarun = load_data('/Volumes/Analysis/2011-06-30-0/d02-08-norefit/data003-from-d02_08/data003-from-d02_08');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun1 = load_data('/Volumes/Analysis/2011-06-30-0/d02-08-norefit/data007-from-d02_08/data007-from-d02_08');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);
datarun1 = load_neurons(datarun1);

% wn_movie_name = 'RGB-2-4-0.48-11111-300x300-60.35.xml'; % for data002
% wn_movie_name = 'BW-2-4-0.48-22222-300x300-60.35.xml'; % for data007
wn_movie_name = 'BW-2-4-0.48-11111-300x300-60.35.xml';
% [inputs, refresh, duration] = get_wn_movie_ath(datarun1, wn_movie_name);


vormap = load('/Volumes/Data/Backup/stimuli/maps/2011-06-30-0_f03_voronoicones/map-0000.txt');
figure
imagesc(vormap)


vorrun = load_data('/Volumes/Analysis/2011-06-30-0/d02-08-norefit/data006-from-d02_08/data006-from-d02_08');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);

% [inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-gaussian-1-4-0.16-11111-1899x1-60.35.xml');

%% find stable cells

tmp_map = vormap;
for i=1:max(tmp_map(:))  
    tmp_map(tmp_map==i) = 0.3+(rand(1)-0.5)/8;
end

cell_type = 4;
bords = 30;
figure
set(gcf, 'position', [58 102 1182 996])
cnt=1;
tmp = zeros(600,600,3);
tmp(:,:,2) = tmp_map;
for i=82:90%find(ww./pixn<5)
    visionID = datarun.cell_types{cell_type}.cell_ids(i);
    cellInd = find(datarun.cell_ids == visionID);
    sta = imresize(double(datarun.stas.stas{cellInd}(:,:,1,4)),600/datarun.stimulus.field_height,'nearest');
    sta1 = imresize(double(datarun1.stas.stas{cellInd}(:,:,1,4)),600/datarun1.stimulus.field_height,'nearest');
 
    
    if datarun.stas.polarities{cellInd} == 1
        [tt, ind] = max(sta(:));
        [tt1, ~] = max(sta1(:));
    else
        [tt, ind] = min(sta(:));
        [tt1, ~] = min(sta1(:));
    end
    [a,b] = ind2sub([600,600],ind);
    
    sta = sta/tt;
    sta1 = sta1/tt1;

    tmp(:,:,1) = sta;
    tmp(:,:,3) = sta1; 
    subplot(3,3,cnt)
    imagesc(tmp)
    axis([b-bords b+bords a-bords a+bords])    
    title(['cell ID ', int2str(cellInd), ', vision ID ', int2str(visionID), ', i=', int2str(i)])
    
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

figure
plot_rf_summaries(datarun, {4}, 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun1, {4}, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')


