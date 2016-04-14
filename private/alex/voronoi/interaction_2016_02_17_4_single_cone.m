%% load stuff

local_path = '/Volumes/Analysis/';


%SC run 1
datarun = load_data([local_path, '2016-02-17-4/d00-05-norefit/data001/data001']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = load_neurons(datarun);

[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-2-6-0.48-11111-400x300-60.35.xml');
inputs = squeeze(inputs);

save('/Users/alexth/Desktop/inputs_tmp.mat', 'inputs', 'refresh', 'duration', '-v7.3')


a = load('/Volumes/Analysis/2016-02-17-4/cone_data/manual/map_data001_manual_postexp_info.mat');
cones = a.cones;
% map = load('/Volumes/Analysis/2016-02-17-4/cone_data/manual/map_data001_manual_postexp.txt');
% figure
% imagesc(map)
% hold on
% plot(cones(:,1), cones(:,2), 'bx');



%% 

% denoise 
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

% cone portraits
cone_table = zeros(length(datarun.cell_ids), length(cones));
for i=1:length(cones)
    a = round(cones(i,:)/2);
    tmp = squeeze(comb_sta(a(2),a(1),:))';
    a=find( abs(tmp)>max_pix*0.4 & abs(tmp)>thresh);
    if ~isempty(a)
        cone_table(a,i) = 1;
    end
end

tmp = raw_sta * pol;

for kkk=  541 

    kkk
%     close all
    datarunID = find(datarun.cell_ids==kkk);
    
    for i=1:length(datarun.cell_types)
        if ~isempty(find(datarun.cell_types{i}.cell_ids==kkk, 1))
            cell_type = datarun.cell_types{i}.name;
            break
        end
    end
    
    
    center_cones = find(cone_table(datarunID,:));
    length(center_cones)
   
    spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/(refresh)); % spikes in frames
    spikes(spikes<100)=[];
    
    
    
    spikes(spikes>size(inputs,2)) = [];
    spikes_tmp = spikes;
    spike_rate=zeros(size(inputs,2),1);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        spikes_tmp(ia)=[];
    end
    clear spikes_tmp
    
    
    
    if length(center_cones)>1 && length(spikes)>1000
        
        cone_coords = cones(center_cones,:);
        
        raw_sta = comb_sta(:,:,datarunID);%*pol(datarunID);

        filt_inputs = zeros(length(center_cones), duration);
        cnt = 1;
        for current_cone=center_cones
            tmp_noise = noise_sta(cone_coords(cnt,2)-1:cone_coords(cnt,2)+1,cone_coords(cnt,1)-1:cone_coords(cnt,1)+1);
          
            tmp_coord = [cone_coords(cnt,1)-1:cone_coords(cnt,1)+1 cone_coords(cnt,2)+k];
            tmp_inputs = cones_inputs{current_cone} - repmat(tmp_noise(:),1, size(cones_inputs{current_cone},2));
            tmp_sta = raw_sta(cone_coords(cnt,2)-1:cone_coords(cnt,2)+1,cone_coords(cnt,1)-1:cone_coords(cnt,1)+1);

            filt_inputs(cnt,:) = sum(tmp_inputs .* repmat(tmp_sta(:),1, size(cones_inputs{current_cone},2)));
            cnt=cnt+1;
        end
        
        
        spikes(spikes>duration) = [];
        spikes_tmp = spikes;
        spike_rate=zeros(duration,1);
        while ~isempty(spikes_tmp)
            [~, ia, ~] = unique(spikes_tmp);
            spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
            spikes_tmp(ia)=[];
        end
        clear spikes_tmp

        err =[];
        try  loglikratio = fit_normal_cdfs(filt_inputs, spikes-1, center_cones');
        catch err
        end
        
        if isempty(err)
            save(['/Volumes/Analysis/2016-02-17-4/cone_data/manual/cell_', int2str(kkk), '.mat'], 'loglikratio');
        end
    end
    cnt1 = cnt1+1;
end


datarunID = find(vorrun.cell_ids==3736);
center_cones = find(cone_table(datarunID,:))'


all_array = zeros((length(center_cones)-1)*4,length(center_cones)*4);
for cone1 = 1:length(center_cones)
    
    for cone2=cone1+1:length(center_cones)
        tmp = loglikratio(cone1, cone2, :);
        tmp=reshape(tmp,3,3);
        all_array(cone1*4-3:cone1*4-1,cone2*4-3:cone2*4-1) = tmp;
    end
end
all_array = all_array/max(abs(all_array(:)))/2+0.5;
all_array = [ones(1,size(all_array,2))-0.5; all_array];
all_array = repmat(all_array,1,1,3);
% all_array(:,:,[1 3]) = 0;
tmp = all_array;
x0 = 0.5;
k = 10;
tmp = ones(size(tmp))./(1+exp(-k*(tmp-x0)));
figure
set(gcf, 'position', [-1822         131        1033         974]);
subplot('position', [0 0 1 1])
imagesc(tmp);
set(gca, 'xtick', 0,'xticklabel','')
set(gca, 'ytick', 0,'yticklabel','')
set(gca, 'visible', 'off')
% set(gca, 'xtick', 2:4:length(center_cones)*4,'xticklabel',int2str(center_cones))
% set(gca, 'ytick', 2:4:length(center_cones)*4,'yticklabel',int2str(center_cones))
xlabel('cone 1')
ylabel('cone 2')
set(gca,'dataaspectratio', [1 1 1])
saveas(gcf, '/Users/alexth/Dropbox/Lab/Transfer/Alex_to_EJ/new_talk/2011-12-13-2/offm_3736/loglik_2d_matrix.svg')
% hold on
% for i = 3.5:3:length(center_cones)*3
%     line([0, length(center_cones)*3+1], [i, i], 'color', [1 1 1]*0.5, 'linewidth',10)
%     line([i, i],[0, length(center_cones)*3+1], 'color', [1 1 1]*0.5, 'linewidth',10)
% end
% for i = 3.5:3:length(center_cones)*3
%     line([0, length(center_cones)*3+1], [i, i], 'color', 'k' )
%     line([i, i],[0, length(center_cones)*3+1], 'color', 'k')
% end
line([0, length(center_cones)*3+1], [0, length(center_cones)*3+1], 'color', 'k')



for kkk= vorrun.cell_types{1}.cell_ids(1:end)  %vorrun.cell_ids(22:end)
    datarunID = find(vorrun.cell_ids==kkk);
    
    for i=1:length(vorrun.cell_types)
        if ~isempty(find(vorrun.cell_types{i}.cell_ids==kkk, 1))
            cell_type = vorrun.cell_types{i}.name;
            break
        end
    end
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID})-noise_vor;
    thresh = mean([robust_std(raw_sta(:,5)),robust_std(raw_sta(:,10))])*5;
    
    a = min(raw_sta(:));
    b = max(raw_sta(:));
    if abs(a)>abs(b) % OFF cell
        pol = -1;
        [~,pos] = find(raw_sta==a,1);
    else % ON cell
        pol = 1;
        [~,pos] = find(raw_sta==b,1);
        a = b;
    end
    tmp = raw_sta * pol;
    
    center_cones = find(sum(tmp(:,pos-1:pos+1)>max(tmp(:))*0.3,2));
    center_cones(isnan(cones(center_cones,1)))=[];
%     center_cones = find(sum(tmp(:,pos-1:pos+1)>thresh,2));
%     center_cones(isnan(cones(center_cones,1)))=[];
    if ~isempty(find(center_cones==441,1)) & ~isempty(find(center_cones==389,1))
        kkk
    end
end

find(vorrun.cell_ids==6991)






center_cones = [390,441]';
cone_cells = find(sum(cone_table(:,center_cones),2)==2);
all_cell_type(cone_cells);
cnt1 = 1;
clear comb_loglikratio
for kkk=  vorrun.cell_ids(cone_cells)
    kkk
    datarunID = find(vorrun.cell_ids==kkk);
    raw_sta = squeeze(vorrun.stas.stas{datarunID})-noise_vor;
    thresh = mean([robust_std(raw_sta(:,5)),robust_std(raw_sta(:,10))])*5;
    
    a = min(raw_sta(:));
    b = max(raw_sta(:));
    if abs(a)>abs(b) % OFF cell
        pol = -1;
        [~,pos] = find(raw_sta==a,1);
    else % ON cell
        pol = 1;
        [~,pos] = find(raw_sta==b,1);
        a = b;
    end
    tmp = raw_sta * pol;
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    
    if abs(a) > thresh*1.5 && length(spikes)>300
        raw_sta = raw_sta(center_cones, end-1:-1:10);
        filt_inputs = zeros(length(center_cones), size(inputs,2)-sta_params.length+1);
        cnt = 1;
        for current_cone=center_cones'
            filt_inputs(cnt,:)=conv(inputs(current_cone,:)-noise_vor(current_cone,1), raw_sta(cnt,:),'valid');
            cnt=cnt+1;
        end
        spikes_tmp = spikes;
        spikes_tmp(spikes<sta_params.length) = [];
        err =[];
        try  loglikratio = fit_normal_cdfs(filt_inputs, spikes_tmp-sta_params.length+1, center_cones);
        catch err
        end
        if isempty(err)
            comb_loglikratio{cnt1} = loglikratio;
        end
        
    end
    cnt1 = cnt1+1;
    
end
save(['/Volumes/Analysis/2011-12-13-2/cone_data/manual/cones_', int2str(center_cones(1)),'_',...
    int2str(center_cones(2)), '.mat'], 'comb_loglikratio');

all_array = [];used_cell = [];
for i = 1:length(comb_loglikratio)
    tmp = squeeze(comb_loglikratio{i});
    if ~isempty(tmp)
        tmp=reshape(tmp(2,:),3,3);
        all_array = [all_array tmp];
        used_cell = [used_cell i];
    end
end

all_array = all_array/max(abs(all_array(:)))/2+0.5;
all_array = repmat(all_array,1,1,3);
% all_array(:,:,[1 3]) = 0;
tmp = all_array;
x0 = 0.5;
k = 10;
tmp = ones(size(tmp))./(1+exp(-k*(tmp-x0)));
figure
set(gcf, 'position', [-1009         722         691         233]);
imagesc(tmp);
set(gca, 'xtick', 2:3:length(used_cell)*3,'xticklabel',int2str(cone_cells(used_cell)))
set(gca, 'ytick', 2,'yticklabel','')
xlabel(['cone ', int2str(center_cones(2))])
ylabel(['cone ', int2str(center_cones(1))])
set(gca,'dataaspectratio', [1 1 1])
hold on
for i = 3.5:3:length(used_cell)*3
    line([i, i],[0, length(used_cell)*3+1], 'color', [1 1 1]*0.5, 'linewidth',10)
end
for i = 3.5:3:length(used_cell)*3
    line([i, i],[0, length(used_cell)*3+1], 'color', 'k')
end









