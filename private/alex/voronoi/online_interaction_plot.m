%% load stuff

local_path = '/Volumes/Acquisition/Analysis/';

% vormap = load(['/Volumes/Analysis/2015-10-29-1/stimuli/maps/map_data001.txt']);

vormap = load(['/Volumes/Analysis/2015-10-29-1/stimuli/maps/map_data004.txt']);

figure
colormap gray
imagesc(vormap)


% vorrun = load_data([local_path, '2015-10-29-1/data003/data003']);
vorrun = load_data([local_path, '2015-10-29-1/data007/data007']);
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_sta(vorrun);
vorrun = load_neurons(vorrun);
tic
[inputs, refresh, duration] = get_wn_movie_ath_rgb(vorrun, 'RGB-1-10-0.48-11111-882x1.xml');
% [inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-882x1.xml');
% [inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-510x1.xml');
toc



%% all on
my_cells = vorrun.cell_types{3}.cell_ids;

for kkk= my_cells
    close all
    datarunID = find(vorrun.cell_ids==kkk);
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
%     figure;
%     plot(raw_sta')
    [full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);
    
    center_cones = find(raw_sta(:,27)>0.02);
    if length(center_cones)>1 && length(center_cones)<25
        
        x = mean(cones(center_cones,1));
        y = mean(cones(center_cones,2));
        tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
        [~, ic] = sort(tmp);
        ic(isnan(tmp(ic))) = [];
        far_cones = ic(end-15:end)';
        
        select_cones = [center_cones; far_cones];
        
        
        voronoi_regions = full_sta(:,:,27);
        voronoi_regions = voronoi_regions/(max(voronoi_regions(:))*1.5);
        comb = voronoi_regions;
        
        sta_params.length = 15;
        sta_params.offset = 0;
        fraction = 0.9;
        
        spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
        spikes(spikes<sta_params.length-sta_params.offset)=[];
        
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(select_cones,1:50000), spikes, fraction, sta_params);
        
        filt_inputs = zeros(length(select_cones), size(inputs,2)-sta_params.length+1);
        cnt = 1;
        for current_cone=select_cones'
            filt_inputs(cnt,:)=conv(inputs(current_cone,:), unbiased_sta(cnt,:),'valid');
            cnt=cnt+1;
        end
        spikes_tmp = spikes;
        spikes_tmp(spikes<sta_params.length) = [];
        
        
        nbins_cone1 = 10;
        nbins_cone2 = 6;
        contrast_response_erm(filt_inputs, spikes_tmp-sta_params.length+1, nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb, datarunID);
    end
end

%% all off
my_cells = vorrun.cell_types{4}.cell_ids;

for kkk= my_cells
    close all
    datarunID = find(vorrun.cell_ids==kkk);
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
%     figure;
%     plot(raw_sta')
    [full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);
    
    center_cones = find(raw_sta(:,27)<-0.1);
    if length(center_cones)>1 && length(center_cones)<20
        
        x = mean(cones(center_cones,1));
        y = mean(cones(center_cones,2));
        tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
        [~, ic] = sort(tmp);
        ic(isnan(tmp(ic))) = [];
        far_cones = ic(end-15:end)';
        
        select_cones = [center_cones; far_cones];
        
        
        voronoi_regions = full_sta(:,:,27);
        voronoi_regions = voronoi_regions/(min(voronoi_regions(:))*1.5);
        comb = voronoi_regions;
        
        sta_params.length = 15;
        sta_params.offset = 0;
        fraction = 0.9;
        
        spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
        spikes(spikes<sta_params.length-sta_params.offset)=[];
        
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(select_cones,1:50000), spikes, fraction, sta_params);
        
        filt_inputs = zeros(length(select_cones), size(inputs,2)-sta_params.length+1);
        cnt = 1;
        for current_cone=select_cones'
            filt_inputs(cnt,:)=conv(inputs(current_cone,:), unbiased_sta(cnt,:),'valid');
            cnt=cnt+1;
        end
        spikes_tmp = spikes;
        spikes_tmp(spikes<sta_params.length) = [];
        
        
        nbins_cone1 = 6;
        nbins_cone2 = 10;
        contrast_response_erm(filt_inputs, spikes_tmp-sta_params.length+1, nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb, datarunID);
    end
end

%% all off parasols
my_cells = vorrun.cell_types{2}.cell_ids;

for kkk= my_cells
    close all
    datarunID = find(vorrun.cell_ids==kkk);
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
%     figure;
%     plot(raw_sta')
    [full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);
    
    center_cones = find(raw_sta(:,27)<-0.03);
    if length(center_cones)>1 && length(center_cones)<20
        
        x = mean(cones(center_cones,1));
        y = mean(cones(center_cones,2));
        tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
        [~, ic] = sort(tmp);
        ic(isnan(tmp(ic))) = [];
        far_cones = ic(end-15:end)';
        
        select_cones = [center_cones; far_cones];
        
        
        voronoi_regions = full_sta(:,:,27);
        voronoi_regions = voronoi_regions/(min(voronoi_regions(:))*1.5);
        comb = voronoi_regions;
        
        sta_params.length = 15;
        sta_params.offset = 0;
        fraction = 0.9;
        
        spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
        spikes(spikes<sta_params.length-sta_params.offset)=[];
        
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(select_cones,1:50000), spikes, fraction, sta_params);
        
        filt_inputs = zeros(length(select_cones), size(inputs,2)-sta_params.length+1);
        cnt = 1;
        for current_cone=select_cones'
            filt_inputs(cnt,:)=conv(inputs(current_cone,:), unbiased_sta(cnt,:),'valid');
            cnt=cnt+1;
        end
        spikes_tmp = spikes;
        spikes_tmp(spikes<sta_params.length) = [];
        
        
        nbins_cone1 = 6;
        nbins_cone2 = 10;
        contrast_response_erm(filt_inputs, spikes_tmp-sta_params.length+1, nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb, datarunID);
    end
end


%% all on parasol
my_cells = vorrun.cell_types{1}.cell_ids;

for kkk= my_cells
    close all
    datarunID = find(vorrun.cell_ids==kkk);
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
%     figure;
%     plot(raw_sta')
    [full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);
    
    center_cones = find(raw_sta(:,27)>0.02);
    if length(center_cones)>1 && length(center_cones)<25
        
        x = mean(cones(center_cones,1));
        y = mean(cones(center_cones,2));
        tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
        [~, ic] = sort(tmp);
        ic(isnan(tmp(ic))) = [];
        far_cones = ic(end-15:end)';
        
        select_cones = [center_cones; far_cones];
        
        
        voronoi_regions = full_sta(:,:,27);
        voronoi_regions = voronoi_regions/(max(voronoi_regions(:))*1.5);
        comb = voronoi_regions;
        
        sta_params.length = 15;
        sta_params.offset = 0;
        fraction = 0.9;
        
        spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
        spikes(spikes<sta_params.length-sta_params.offset)=[];
        
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(select_cones,1:50000), spikes, fraction, sta_params);
        
        filt_inputs = zeros(length(select_cones), size(inputs,2)-sta_params.length+1);
        cnt = 1;
        for current_cone=select_cones'
            filt_inputs(cnt,:)=conv(inputs(current_cone,:), unbiased_sta(cnt,:),'valid');
            cnt=cnt+1;
        end
        spikes_tmp = spikes;
        spikes_tmp(spikes<sta_params.length) = [];
        
        
        nbins_cone1 = 10;
        nbins_cone2 = 6;
        contrast_response_erm(filt_inputs, spikes_tmp-sta_params.length+1, nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb, datarunID);
    end
end


%% 
my_cells = vorrun.cell_types{1}.cell_ids;

max_cone = max(vormap(:));
for kkk= my_cells
%     close all
    datarunID = find(vorrun.cell_ids==kkk);
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
  full_sta=zeros(size(vormap,1),size(vormap,2),size(raw_sta,2));

    for i=1:max_cone
        [a, b] = find(vormap==i);
        
        if ~isempty(a)
            cones(i,1) = mean(a);
            cones(i,2) = mean(b);
            for j = 1:length(a)
                full_sta(a(j),b(j),:) = raw_sta(i,:,30);
            end
        else
            cones(i,:) = nan;
        end
    end
    tmp = full_sta/max(abs(full_sta(:)))/2+0.5;
    figure
    imagesc(tmp)
    drawnow
end