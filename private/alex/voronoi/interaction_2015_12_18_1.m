%% load stuff

% map
vormap = load('/Volumes/Analysis/2015-12-18-1/stimuli/maps/map_data002.txt');
figure
colormap gray
imagesc(vormap)

temp = uint8(zeros(size(vormap)));
voronoi_contours = cell(max(vormap(:)),2);
figure
colormap gray
td = vormap;
td(vormap>0) = 1;
imagesc(td)
hold on
for i=1:max(vormap(:))
    tmpmap=temp;
    if ~isempty(find(vormap==i,1))
        tmpmap(vormap==i)=1;
        dd = imresize(tmpmap,5,'method', 'nearest');
        [r c] = find(dd,1);
        contour = bwtraceboundary(dd,[r c],'W',8,Inf,'counterclockwise');
        contour= round(contour/5)+0.5;
        plot(contour(:,2), contour(:,1), 'linewidth', 2)
        voronoi_contours{i, 1} = contour(:,2);
        voronoi_contours{i, 2} = contour(:,1);
    end
end

%SC run 1
datarun = load_data('/Volumes/Analysis/2015-12-18-1/d00-12-norefit/data002/data002');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);

%SC run 2
datarun2 = load_data('/Volumes/Analysis/2015-12-18-1/d00-12-norefit/data005/data005');
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2);

% voronoi run
vorrun = load_data('/Volumes/Analysis/2015-12-18-1/d00-12-norefit/data004/data004');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_sta(vorrun);
vorrun = load_neurons(vorrun);
[inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-569x1-119.5.xml');


%% plot stuff
save_path = '/Volumes/Analysis/2015-12-18-1/data004_d00-12/';
nbins_cone1 = 6;
nbins_cone2 = 6;
sta_params.length = 20;
sta_params.offset = 0;
fraction = 0.9;



for kkk= vorrun.cell_ids(22:end)
    close all
    datarunID = find(vorrun.cell_ids==kkk);
    
    for i=1:length(vorrun.cell_types)
        if ~isempty(find(vorrun.cell_types{i}.cell_ids==kkk))
            cell_type = vorrun.cell_types{i}.name;
            break
        end
    end
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
    %     figure;
    %     plot(raw_sta')
    thresh = mean([robust_std(raw_sta(:,5)),robust_std(raw_sta(:,10))])*5;
    
    a = min(raw_sta(:));
    b = max(raw_sta(:));
    if abs(a)>abs(b) % OFF cell
        pol = -1;
        [~,pos] = find(raw_sta==a,1);
    else % ON cell
        pol = 1;
    end
    tmp = raw_sta * pol;
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    
    if abs(a) > thresh*1.5 && length(spikes)>300
        center_cones = find(sum(tmp(:,pos-1:pos+1)>thresh,2));
        center_cones(isnan(cones(center_cones,1)))=[];
        
        
        %     figure
        %     plot(raw_sta(center_cones,:)')
        timecourse = mean(tmp(center_cones,10:end));
        %     figure
        %     plot(timecourse)
        
        conv_sta = sum(tmp(:,10:end).*repmat(timecourse,size(tmp,1),1),2);
        if length(center_cones)>30
            [~, a] = sort(conv_sta);
            center_cones = sort(a(end-29:end));
        end
        timecourse = mean(tmp(center_cones,10:end));
        conv_sta = sum(tmp(:,10:end).*repmat(timecourse,size(tmp,1),1),2);
        
        % prepare voronoi and single cone STAs
        [full_sta, cones] = expand_voronoi_sta(conv_sta, vormap);
        full_sta = pol*full_sta/(max(full_sta(:))*1.5);
        sta1 = squeeze(datarun.stas.stas{datarunID});
        sta1 = sta1(:,:,5);
        sta2 = squeeze(datarun2.stas.stas{datarunID});
        sta2 = sta2(:,:,5);
        
        if length(center_cones)>1
            
           
            [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(center_cones,1:50000), spikes, fraction, sta_params);
            
            %         figure
            %         plot(unbiased_sta')
            
            filt_inputs = zeros(length(center_cones), size(inputs,2)-sta_params.length+1);
            cnt = 1;
            for current_cone=center_cones'
                filt_inputs(cnt,:)=conv(inputs(current_cone,:), unbiased_sta(cnt,:),'valid');
                cnt=cnt+1;
            end
            spikes_tmp = spikes;
            spikes_tmp(spikes<sta_params.length) = [];
            
            offline_contrast_response(filt_inputs, spikes_tmp-sta_params.length+1,...
                nbins_cone1, nbins_cone2, center_cones, vormap, cones, full_sta, datarunID,...
                kkk, [save_path int2str(nbins_cone2) '_bins/', cell_type '/'], sta1, sta2, voronoi_contours, pol);
        end
    end
end

