datarun = load_data('/Volumes/Analysis/2011-05-11-6/d03-09-norefit/data003-from-d03_09/data003-from-d03_09');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);


datarun1 = load_data('/Volumes/Analysis/2011-05-11-6/d03-09-norefit/data009-from-d03_09/data009-from-d03_09');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);

[inputs, refresh, duration] = get_wn_movie_ath_rgb(datarun, 'RGB-2-4-0.5-11111-300x300-60.35.xml');

vorrun = load_data('/Volumes/Analysis/2011-05-11-6/d03-09-norefit/data005-from-d03_09/data005-from-d03_09'); % RGB
%vorrun = load_data('/Volumes/Analysis/2011-05-11-6/d03-09-norefit/data006-from-d03_09/data006-from-d03_09'); % BW
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);

[inputs_v, refresh_v, duration_v] = get_wn_movie_ath_rgb(vorrun, 'RGB-1-1-0.50-11111-2027x1-60.35.xml'); % RGB
% [inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-1-4-0.50-11111-2027x1-60.35.xml'); % BW

vormap = load('/Volumes/Data/2011-05-11-6/Visual/voronoi_300x300_radius8_data002/map-0000.txt');
figure
imagesc(vormap)
max_ind = max(vormap(:));

tmp_map = vormap;
for i=1:max(tmp_map(:))  
    tmp_map(tmp_map==i) = 0.3+(rand(1)-0.5)/8;
end

offm = datarun.cell_types{4}.cell_ids;
offp = datarun.cell_types{2}.cell_ids;
sbc = datarun.cell_types{5}.cell_ids;
bords = 30;
figure
set(gcf, 'position', [58 102 1182 996])
cnt=1;
tmp = zeros(300,300,3);
tmp(:,:,2) = tmp_map;
for i=1:12
    cellInd = find(datarun.cell_ids == offm(i));
    sta = double(datarun.stas.stas{cellInd}(:,:,2,4));
    sta1 = double(datarun1.stas.stas{cellInd}(:,:,1,4));
    
    if datarun.stas.polarities{cellInd} == 1
        [tt, ind] = max(sta(:));
        [tt1, ~] = max(sta1(:));
    else
        [tt, ind] = min(sta(:));
        [tt1, ~] = min(sta1(:));
    end
    [a,b] = ind2sub([300,300],ind);
    
    sta = sta/tt;
%     sta(sta>0.25) = 1;

    sta1 = sta1/tt1;
%     sta1(sta1>0.25) = 1;

    tmp(:,:,1) = sta;
    tmp(:,:,3) = sta1;    
    subplot(3,4,cnt)
    imagesc(tmp)
    axis([b-bords b+bords a-bords a+bords])
    title(['cell ', int2str(cellInd), ' vision ID ', int2str(offm(i))])
    cnt = cnt+1;
end


clear bords sta sta1 ind tt tt1 cellInd a b

offset = 0;
sta_length = 15;
my_sta=zeros(max_ind,3, sta_length, length(cells));
spike_rate=zeros(duration_v,length(cells)); 
cnt = 1;
for i = cells
    cellInd = find(datarun.cell_ids == offm(i));
    spikes=ceil((vorrun.spikes{cellInd}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spikes(spikes<sta_length-offset)=[];
    
    while ~isempty(spikes)
        [~, ia, ~] = unique(spikes);
        spike_rate(spikes(ia), cnt)=spike_rate(spikes(ia), cnt)+1;
        for j=1:sta_length
            my_sta(:,:,sta_length-j+1,cnt)=my_sta(:,:,sta_length-j+1,cnt)+sum(inputs_v(:,:,spikes(ia)-sta_length+j+offset),3);
        end
        spikes(ia)=[];
    end
    my_sta(:,:,:,cnt)=my_sta(:,:,:,cnt)/sum(spike_rate(:,cnt));
    cnt = cnt+1;
end

tmp_map = vormap;

my_cell = 1;
tt=0:max_ind;
voronoi_sta=zeros(300,300,3,sta_length);
for i=1:max_ind
    [a, b] = find(vormap==tt(i+1));
    if ~isempty(a)
        for j = 1:length(a)
            voronoi_sta(a(j),b(j),:,:) = my_sta(i,:,:, my_cell);
        end
    end
end

tmp = voronoi_sta(:);
abs_max =  0.5/max(abs(tmp));
min(tmp)
max(tmp)

figure
for i=1:9
    subplot(3,3,10-i)
%     colormap gray
    imagesc(voronoi_sta(:,:,:,i)*abs_max+0.5)
    title(['frame back ', int2str(i)])
%     axis([300 365 460 520])
end


%% Plot stuff

% first do Voronoi - simply 3rd frame, no temporal summation
sta_params.length = 7;
sta_params.offset = 0;
offset = 0;
sta_length=4;
filepath = '/Users/alexth/Desktop/voronoi/2011-05-11-6/d03-09-norefit/';

for datarunID = 1:length(datarun.cell_ids)
    visionID = datarun.cell_ids(datarunID);
    [folder, my_type] = find_cell_type(datarun, visionID);
    
    if ~exist([filepath,folder],'dir')
        mkdir([filepath,folder]);
    end
    
    
    fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[1 24 1886 1074]);
    
    my_sta=zeros(size(inputs_v,1),3,sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    
    % 108621 frames is 1800s;
    k = [sta_params.length-sta_params.offset, 36207, 72414, 108621];

    for cnt1 = 1:3
        spikes_tmp = spikes(spikes>k(cnt1) & spikes <=k(cnt1+1));
        
        while ~isempty(spikes_tmp)
            [~, ia, ~] = unique(spikes_tmp);
            for j=1:sta_params.length
                my_sta(:,:,sta_params.length-j+1) = my_sta(:, :, sta_params.length-j+1) +...
                    sum(inputs_v(:,:,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),3);
            end
            spikes_tmp(ia)=[];
        end
        % my_sta=my_sta/nspikes;
        
        % get voronoi map sta
        tmp_map = vormap;
        tt=0:max(tmp_map(:));
        vorsta=zeros(300,300,3,sta_params.length);
        coneX = zeros(max(tmp_map(:)),1);
        coneY = coneX;
        cnt = 1;
        for i=1:2027
            [a, b] = find(tmp_map==tt(i+1));
            % find center of this cone
            coneX(i) = mean(a);
            coneY(i) = mean(b);
            if ~isempty(a)
                for j = 1:length(a)
                    vorsta(a(j),b(j),:,:) = my_sta(cnt,:,:);
                end
            end
            cnt=cnt+1;
        end
        vorsta_tmp = vorsta(:, :, :, 5);
        
        tmp = vorsta_tmp(:);
        abs_max =  0.5/max(abs(tmp));

        
        subplot(2,3,cnt1)
        colormap gray
        imagesc(vorsta_tmp*abs_max+0.5)
        title([int2str((k(cnt1+1))/(1000/refresh_v)), ' s'])
        drawnow
    end
    
    
    subplot(2,3,2)
    title(['Vision ID ', int2str(visionID), ', datarun ID ', int2str(datarunID), ', ', folder])

    
    % do single cone WN
    
    sta = double(datarun.stas.stas{datarunID});
    
    tmp = sta(:);
    abs_max =  0.5/max(abs(tmp));
    
    subplot(2,3,6)
    imagesc(sta(:,:,:,4)*abs_max+0.5)
    drawnow    
    
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(datarunID)]));
    close(fig)

end


