%% load stuff

local_path = '/Volumes/Analysis/';


%SC run 1
datarun = load_data([local_path, '2016-02-17-4/d00-05-norefit/data001/data001']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = load_neurons(datarun);

a = load('/Volumes/Analysis/2016-02-17-4/cone_data/manual/map_data001_manual_postexp1_info.mat');
cones = round(a.cones/2);

if 0
    tic
    [full_inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-2-6-0.48-11111-400x300-60.35.xml');
    toc
    
%     full_inputs = reshape(full_inputs, 120000, duration);
%     
%     for i=1:1000:120000
%         inputs = full_inputs(i:i+999,:);
%         save(['/Volumes/Analysis/2016-02-17-4/subunits/inputs_parts/part_',int2str(i),'.mat'], 'inputs', 'refresh', 'duration')
%     end
%     
%     a = load('/Volumes/Analysis/2016-02-17-4/cone_data/manual/map_data001_manual_postexp1_info.mat');
%     cones = round(a.cones/2);
    
    
    full_inputs = reshape(full_inputs, 300,400, duration);
    
    cones_inputs = cell(1,length(cones));
    for i = 1:length(cones)
        cones_inputs{i} = full_inputs(cones(i,2)-1:cones(i,2)+1,cones(i,1)-1:cones(i,1)+1,:);
    end
    
    for i=1:10:length(cones)
        if i+9>length(cones)
            inputs = cones_inputs(i:end);
        else
            inputs = cones_inputs(i:i+9);
        end
        save(['/Volumes/Analysis/2016-02-17-4/subunits/indiv_cone_inputs/cone_',int2str(i),'.mat'], 'inputs', 'refresh', 'duration')
    end
    
    %     a = cones_inputs{474};
    %     b = inputs{4};
    %     tmp = b(:,:,1)';
    %
end
%%

tmp = 0;
for i=1:length(datarun.cell_ids)
    tmp = tmp + squeeze(datarun.stas.stas{i}(:,:,:,4));
end
noise_sta = tmp/length(datarun.cell_ids);

% comb STA
cnt = 1;
clear comb_sta max_pix thresh pol
for kkk=datarun.cell_ids
    datarunID = find(datarun.cell_ids==kkk);
    
    for i=1:length(datarun.cell_types)
        if ~isempty(find(datarun.cell_types{i}.cell_ids==kkk, 1))
            all_cell_type(cnt) = i;
            break
        end
    end
    raw_sta = squeeze(datarun.stas.stas{datarunID}(:,:,:,4))-noise_sta;
    thresh(cnt) = mean([robust_std(raw_sta(:,1)),robust_std(raw_sta(:,6))])*5;
    max_pix(cnt) = max([max(raw_sta(:)) abs(min(raw_sta(:)))]);
    comb_sta(:,:,cnt) = raw_sta;
    
    a = min(raw_sta(:));
    b = max(raw_sta(:));
    if abs(min(raw_sta(:)))>max(raw_sta(:)) % OFF cell
        pol(cnt) = -1;
    else % ON cell
        pol(cnt) = 1;
    end
    
    cnt=cnt+1;
end


% cone connectivity table
cone_table = zeros(length(datarun.cell_ids), length(cones));
for i=1:length(cones)
    tmp = squeeze(comb_sta(cones(i,2),cones(i,1),:))';
    a=find( abs(tmp)>max_pix*0.2 & abs(tmp)>thresh);
    if ~isempty(a)
        cone_table(a,i) = 1;
    end
end


% for kkk =  [1921 3031 3721 4022 4353 4396 4486 4876 5056]
%     datarunID = find(datarun.cell_ids==kkk);
%      center_cones = find(cone_table(datarunID,:));
%     length(center_cones)
% end

% [541 571 586 783 1921 2012 2091 2149 3031 3274 3721 4022 4353 4396 4486 4489 4774 4876 4924 4925  5056 5898 6889]
for kkk =  4876
    
    kkk
    
    datarunID = find(datarun.cell_ids==kkk);
    
    center_cones = find(cone_table(datarunID,:));
    length(center_cones)
    
    % get duration and refresh
    load('/Volumes/Analysis/2016-02-17-4/subunits/indiv_cone_inputs/cone_1.mat')
    
    spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/(refresh)); % spikes in frames
    spikes(spikes<10)=[];
    spikes(spikes>duration) = [];
    spikes = spikes-1;
    
    spikes_tmp = spikes;
    spike_rate=zeros(duration,1);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        spikes_tmp(ia)=[];
    end
    clear spikes_tmp
    
    if length(center_cones)>1 && length(spikes)>1000
        
        cone_coords = cones(center_cones,:);
        center_cones_inputs = cell(1,length(center_cones));
        % load inputs
        cnt = 1;
        for i = center_cones
            if mod(i,10)==0
                tt = floor(i/10)*10+1-10;
            else
                tt = floor(i/10)*10+1;
            end
            p = int2str(tt);
            load(['/Volumes/Analysis/2016-02-17-4/subunits/indiv_cone_inputs/cone_',p,'.mat'])
            center_cones_inputs{cnt} = reshape(inputs{i-tt+1}, 9, duration);
            cnt=cnt+1;
        end
        
        raw_sta = comb_sta(:,:,datarunID);
        
        filt_inputs = zeros(length(center_cones), duration);
        cnt = 1;
        for current_cone=center_cones
            tmp_noise = noise_sta(cone_coords(cnt,2)-1:cone_coords(cnt,2)+1,cone_coords(cnt,1)-1:cone_coords(cnt,1)+1);
            tmp_inputs = center_cones_inputs{cnt} - repmat(tmp_noise(:),1, size(center_cones_inputs{cnt},2));
            tmp_sta = raw_sta(cone_coords(cnt,2)-1:cone_coords(cnt,2)+1,cone_coords(cnt,1)-1:cone_coords(cnt,1)+1);
            
            filt_inputs(cnt,:) = sum(tmp_inputs .* repmat(tmp_sta(:),1, size(center_cones_inputs{cnt},2)));
            cnt=cnt+1;
        end
        
%         path2save = ['/Volumes/Analysis/2016-02-17-4/subunits/surface_plots/cell_', int2str(kkk),'/'];
%         if ~isdir(path2save)
%             mkdir(path2save);
%         end
%         
%         err1 = [];
%         try load(['/Volumes/Analysis/2016-02-17-4/cone_data/manual/cell_',int2str(kkk),'.mat'])
%         catch err1
%         end
%         LLR = loglikratio;
       
        err =[];
        try  [loglikratio, res] = fit_normal_cdfs(filt_inputs, spike_rate, center_cones');%, kkk, LLR);
        catch err
        end
        if isempty(err)
            save(['/Volumes/Analysis/2016-02-17-4/cone_data/manual/logistic/cell_', int2str(kkk), '.mat'], 'loglikratio', 'res', 'center_cones', 'raw_sta');
        end
        

        if 0
            figure
            colormap gray
            imagesc(raw_sta)
            hold on
            for i = center_cones
                text(cones(i,1), cones(i,2), int2str(i), 'color', 'r', 'fontsize', 12)
            end
            axis([min(cones(center_cones,1))-10 max(cones(center_cones,1))+10 ...
                min(cones(center_cones,2))-10 max(cones(center_cones,2))+10]);
            saveas(gcf, [path2save, 'raw_sta.bmp']);
            close all
            
            
            if isempty(err1)
                llr_subset = 1:length(center_cones);
                all_array = zeros((length(center_cones)-1)*6,length(center_cones)*6);
                
                for cone1 = 1:length(center_cones)
                    
                    for cone2=cone1+1:length(center_cones)
                        tmp = LLR(llr_subset(cone1), llr_subset(cone2), :);
                        tmp=reshape(tmp,5,5);
                        all_array(cone1*6-5:cone1*6-1,cone2*6-5:cone2*6-1) = tmp;
                    end
                end
                all_array = all_array/max(abs(all_array(:)))/2+0.5;
                all_array = [ones(1,size(all_array,2))-0.5; all_array];
                all_array = repmat(all_array,1,1,3);
                tmp = all_array;
                x0 = 0.5;
                k = 10;
                tmp = ones(size(tmp))./(1+exp(-k*(tmp-x0)));
                figure
                set(gcf, 'position', [-1822         131        1033         974]);
                imagesc(tmp(:,6:end,:));
                set(gca, 'xtick', 0,'xticklabel','')
                set(gca, 'ytick', 0,'yticklabel','')
                set(gca, 'xtick', 4:6:length(center_cones(2:end))*6,'xticklabel',int2str(center_cones(2:end)'), 'fontsize', 10)
                set(gca, 'ytick', 4:6:length(center_cones)*6,'yticklabel',int2str(center_cones'))
                xlabel('cone 1')
                ylabel('cone 2')
                set(gca,'dataaspectratio', [1 1 1])
                title('cdf, resampling with duplication')
                saveas(gcf, [path2save, 'LLR.bmp']);
                close all
            end
        end

    end
end


