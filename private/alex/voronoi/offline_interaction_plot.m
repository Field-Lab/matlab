%% load stuff

local_path = '/Volumes/Analysis/';

vormap = load(['/Volumes/Analysis/2015-10-29-1/stimuli/maps/map_data001.txt']);



figure
colormap gray
imagesc(vormap)


% vorrun = load_data([local_path, '2015-10-29-1/data003/data003']);
% vorrun = load_data([local_path, '2015-11-09-4/data006/data006']);
% vorrun = load_data([local_path, '2016-01-05-1/streamed/data003/data003']);
vorrun = load_data([local_path, '2015-12-18-1/d00-12-norefit/data004/data004']);
vorrun = load_data([local_path, '2015-11-09-4/d00-8-norefit/data006/data006']);

vorrun = load_data([local_path, '2015-10-29-1/d00-12-norefit/data005/data005']);
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_sta(vorrun);
vorrun = load_neurons(vorrun);
tic
[inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-882x1.xml');
% [inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-772x1-119.5.xml');
% [inputs, refresh, duration] = get_wn_movie_ath_rgb(vorrun, 'RGB-1-10-0.48-11111-882x1.xml');
% [inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-882x1.xml');
% [inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-510x1.xml');
toc

%% manual map
save_path = '/Volumes/Analysis/2015-10-29-1/data005_d00-12/ofp_correct_';

my_cells = vorrun.cell_types{2}.cell_ids;
pol = -1;
nbins_cone1 = 6;
nbins_cone2 = 6;

sta1=25;
sta2=26;

for kkk= my_cells(1:end)
    close all
    datarunID = find(vorrun.cell_ids==kkk);
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
%     figure;
%     plot(raw_sta')
    [full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);
    
    thresh = mean([robust_std(raw_sta(:,5)),robust_std(raw_sta(:,10))])*4;

    if pol==1
        center_cones = find(raw_sta(:,sta1)>thresh | raw_sta(:,sta2)>thresh)-1; %%%%%%%%%
        if sum(raw_sta(:,sta1)>thresh)>sum(raw_sta(:,sta2)>thresh)
            sta_frame = sta1;
        else
            sta_frame = sta2;
        end
    else
        center_cones = find(raw_sta(:,sta1)<-thresh | raw_sta(:,sta2)<-thresh)-1; %%%%%%%%%
        if sum(raw_sta(:,sta1)<-thresh)>sum(raw_sta(:,sta2)<-thresh)
            sta_frame = sta1;
        else
            sta_frame = sta2;
        end
    end
    center_cones(isnan(cones(center_cones,1)))=[];
    if length(center_cones)>1 && length(center_cones)<20
        
        x = mean(cones(center_cones,1));
        y = mean(cones(center_cones,2));
        tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
        [~, ic] = sort(tmp);
        ic(isnan(tmp(ic))) = [];
        far_cones = ic(end-15:end)';
        far_cones(far_cones==length(cones))=[];
        far_cones(isnan(cones(far_cones,1)))=[];
        select_cones = [center_cones; far_cones];
        
        voronoi_regions = full_sta(:,:,sta_frame);
        voronoi_regions = voronoi_regions/(min(voronoi_regions(:))*1.5);
        comb = voronoi_regions;
        
        sta_params.length = 15;
        sta_params.offset = 0;
        fraction = 0.9;
        
        spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
        spikes(spikes<sta_params.length-sta_params.offset)=[];
        
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(select_cones+1,1:100000), spikes, fraction, sta_params);

        filt_inputs = zeros(length(select_cones), size(inputs,2)-sta_params.length+1);
        cnt = 1;
        for current_cone=select_cones'
            filt_inputs(cnt,:)=conv(inputs(current_cone+1,:), unbiased_sta(cnt,:),'valid');
            cnt=cnt+1;
        end
        spikes_tmp = spikes;
        spikes_tmp(spikes<sta_params.length) = [];
        
        online_contrast_response(filt_inputs, spikes_tmp-sta_params.length+1,...
            nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb, datarunID,...
            kkk, [save_path int2str(nbins_cone2) '_bins/']);
    end
end



%% bayesian map

save_path = '/Volumes/Analysis/2015-10-29-1/data005_d00-12/ofp_correct_';

my_cells = vorrun.cell_types{2}.cell_ids;
pol = -1;
nbins_cone1 = 6;
nbins_cone2 = 6;

sta1=25;
sta2=26;

for kkk= my_cells(1:end)
    close all
    datarunID = find(vorrun.cell_ids==kkk);
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
%     figure;
%     plot(raw_sta')
    [full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);
    
    thresh = mean([robust_std(raw_sta(:,5)),robust_std(raw_sta(:,10))])*4;

    if pol==1
        center_cones = find(raw_sta(:,sta1)>thresh | raw_sta(:,sta2)>thresh)-1; %%%%%%%%%
        if sum(raw_sta(:,sta1)>thresh)>sum(raw_sta(:,sta2)>thresh)
            sta_frame = sta1;
        else
            sta_frame = sta2;
        end
    else
        center_cones = find(raw_sta(:,sta1)<-thresh | raw_sta(:,sta2)<-thresh); %%%%%%%%%
        if sum(raw_sta(:,sta1)<-thresh)>sum(raw_sta(:,sta2)<-thresh)
            sta_frame = sta1;
        else
            sta_frame = sta2;
        end
    end
    center_cones(isnan(cones(center_cones,1)))=[];
    if length(center_cones)>1 && length(center_cones)<20
        
        x = mean(cones(center_cones,1));
        y = mean(cones(center_cones,2));
        tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
        [~, ic] = sort(tmp);
        ic(isnan(tmp(ic))) = [];
        far_cones = ic(end-15:end)';
        far_cones(far_cones==length(cones))=[];
        far_cones(isnan(cones(far_cones,1)))=[];
        select_cones = [center_cones; far_cones];
        
        voronoi_regions = full_sta(:,:,sta_frame);
        voronoi_regions = voronoi_regions/(min(voronoi_regions(:))*1.5);
        comb = voronoi_regions;
        
        sta_params.length = 15;
        sta_params.offset = 0;
        fraction = 0.9;
        
        spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
        spikes(spikes<sta_params.length-sta_params.offset)=[];
        
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(select_cones,1:100000), spikes, fraction, sta_params);

        filt_inputs = zeros(length(select_cones), size(inputs,2)-sta_params.length+1);
        cnt = 1;
        for current_cone=select_cones'
            filt_inputs(cnt,:)=conv(inputs(current_cone,:), unbiased_sta(cnt,:),'valid');
            cnt=cnt+1;
        end
        spikes_tmp = spikes;
        spikes_tmp(spikes<sta_params.length) = [];
        
        online_contrast_response(filt_inputs, spikes_tmp-sta_params.length+1,...
            nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb, datarunID,...
            kkk, [save_path int2str(nbins_cone2) '_bins/']);
    end
end


