a = rand(100)>0.5;

figure
colormap gray

set(gca,'position', [0 0 1 1])
imagesc(a)

set(gca,'dataaspectratio', [1 1 1])
axis off


vormap = load([local_path, '2011-12-13-2/Visual/2011-12-13-2_f04_vorcones/map-0000.txt']);

vormap(vormap>0) = 1;
figure
colormap gray
set(gca,'position', [0 0 1 1])
imagesc(vormap)
set(gca,'dataaspectratio', [1 1 1])
axis off


sta = squeeze(datarun.stas.stas{168});
sta_snippet = imresize(double(sta(:,:,4)), 2, 'nearest');


colormap gray
set(gca,'position', [0 0 1 1])
imagesc(sta_snippet)
set(gca,'dataaspectratio', [1 1 1])
axis off
axis([290 355 470 515]+1)

%% coarse sta

datarun = load_data('/Volumes/My Passport for Mac/work/data/2015-03-09-2/d05-27-norefit/data014-from-d05-d27/data014-from-d05-d27');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');

sta = squeeze(datarun.stas.stas{23});
figure
colormap gray
imagesc(sta(:,:,27))
set(gca,'position', [0 0 1 1])
set(gca,'dataaspectratio', [1 1 1])
axis off
axis([22 29 19 24])

%% combination of voronoi and cones

comb2 = comb;
comb1 = comb(:,:,1);
comb1(comb1<0.3) = 0;
comb2(:,:,1) = comb1*1.5;
comb1 = comb(:,:,3);
comb1(comb1<0.3) = 0;
comb2(:,:,3) = comb1*1.5;
comb1 = comb(:,:,2);
comb1(comb1<0.3) = 0;
comb2(:,:,2) = comb1;
figure
imagesc(comb2)
set(gca,'position', [0 0 1 1])
set(gca,'dataaspectratio', [1 1 1])
axis off
axis([290 355 470 515]+1)


%% ID 168
datarunID = 168;

visionID = datarun.cell_ids(datarunID);
raw_sta = squeeze(vorrun.stas.stas{datarunID});
[full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);

center_cones = find(raw_sta(:,27)<-0.07);

x = mean(cones(center_cones,1));
y = mean(cones(center_cones,2));
tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
[~, ic] = sort(tmp);
ic(isnan(tmp(ic))) = [];
far_cones = ic(end-15:end)';

select_cones = [center_cones; far_cones];


sta = squeeze(datarun.stas.stas{datarunID});
sta1 = squeeze(datarun1.stas.stas{datarunID});
sta_snippet = -imresize(double(sta(:,:,4)), 2, 'nearest');
sta1_snippet = -imresize(double(sta1(:,:,4)), 2, 'nearest');
voronoi_regions = full_sta(:,:,27);
voronoi_regions = voronoi_regions/(min(voronoi_regions(:))*1.5);
comb = zeros(600,600,3);
comb(:,:,1) = sta_snippet/max(sta_snippet(:));
comb(:,:,3) = sta1_snippet/max(sta1_snippet(:));
comb(:,:,2) = voronoi_regions;

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


nbins_cone1 = 4;
nbins_cone2 = 4;
contrast_response_erm(filt_inputs, spikes_tmp-sta_params.length+1, nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb, datarunID);


%% ID 271
datarunID = 271;

visionID = datarun.cell_ids(datarunID);
raw_sta = squeeze(vorrun.stas.stas{datarunID});
[full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);

center_cones = find(raw_sta(:,27)>0.05);

x = mean(cones(center_cones,1));
y = mean(cones(center_cones,2));
tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
[~, ic] = sort(tmp);
ic(isnan(tmp(ic))) = [];
far_cones = ic(end-15:end)';

select_cones = [center_cones; far_cones];


sta = squeeze(datarun.stas.stas{datarunID});
sta1 = squeeze(datarun1.stas.stas{datarunID});
sta_snippet = imresize(double(sta(:,:,4)), 2, 'nearest');
sta1_snippet = imresize(double(sta1(:,:,4)), 2, 'nearest');
voronoi_regions = full_sta(:,:,27);
voronoi_regions = voronoi_regions/(max(voronoi_regions(:))*1.5);
comb = zeros(600,600,3);
comb(:,:,1) = sta_snippet/max(sta_snippet(:));
comb(:,:,3) = sta1_snippet/max(sta1_snippet(:));
comb(:,:,2) = voronoi_regions;

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

nbins_cone1 = 4;
nbins_cone2 = 4;
contrast_response_erm(filt_inputs, spikes_tmp-sta_params.length+1, nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb, datarunID);



%% all off
my_cells = datarun.cell_types{4}.cell_ids;

for kkk= my_cells(62:end)
    close all
    datarunID = find(datarun.cell_ids==kkk);
    
    visionID = datarun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
    [full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);
    
    center_cones = find(raw_sta(:,27)<-0.05);
    if length(center_cones)>1 && length(center_cones)<15
        
        x = mean(cones(center_cones,1));
        y = mean(cones(center_cones,2));
        tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
        [~, ic] = sort(tmp);
        ic(isnan(tmp(ic))) = [];
        far_cones = ic(end-15:end)';
        
        select_cones = [center_cones; far_cones];
        
        
        sta = squeeze(datarun.stas.stas{datarunID});
        sta1 = squeeze(datarun1.stas.stas{datarunID});
        sta_snippet = -imresize(double(sta(:,:,4)), 2, 'nearest');
        sta1_snippet = -imresize(double(sta1(:,:,4)), 2, 'nearest');
        voronoi_regions = full_sta(:,:,27);
        voronoi_regions = voronoi_regions/(min(voronoi_regions(:))*1.5);
        comb = zeros(600,600,3);
        comb(:,:,1) = sta_snippet/max(sta_snippet(:));
        comb(:,:,3) = sta1_snippet/max(sta1_snippet(:));
        comb(:,:,2) = voronoi_regions;
        
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
        
        
        nbins_cone1 = 4;
        nbins_cone2 = 4;
        contrast_response_erm(filt_inputs, spikes_tmp-sta_params.length+1, nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb, datarunID);
    end
end


%% all on
my_cells = datarun.cell_types{3}.cell_ids;

for kkk= my_cells(31:end)
    close all
    datarunID = find(datarun.cell_ids==kkk);
    
    visionID = datarun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
    [full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);
    
    center_cones = find(raw_sta(:,27)>0.05);
    if length(center_cones)>1 && length(center_cones)<15
        
        x = mean(cones(center_cones,1));
        y = mean(cones(center_cones,2));
        tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
        [~, ic] = sort(tmp);
        ic(isnan(tmp(ic))) = [];
        far_cones = ic(end-15:end)';
        
        select_cones = [center_cones; far_cones];
        
        
        sta = squeeze(datarun.stas.stas{datarunID});
        sta1 = squeeze(datarun1.stas.stas{datarunID});
        sta_snippet = imresize(double(sta(:,:,4)), 2, 'nearest');
        sta1_snippet = imresize(double(sta1(:,:,4)), 2, 'nearest');
        voronoi_regions = full_sta(:,:,27);
        voronoi_regions = voronoi_regions/(max(voronoi_regions(:))*1.5);
        comb = zeros(600,600,3);
        comb(:,:,1) = sta_snippet/max(sta_snippet(:));
        comb(:,:,3) = sta1_snippet/max(sta1_snippet(:));
        comb(:,:,2) = voronoi_regions;
        
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
        
        
        nbins_cone1 = 4;
        nbins_cone2 = 4;
        contrast_response_erm(filt_inputs, spikes_tmp-sta_params.length+1, nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb, datarunID);
    end
end

%% model
ncones = 10;
sample_size = 10000;
cones = zeros(ncones, sample_size);
for i=1:ncones
    cones(i, :) = rand(1, sample_size)-0.5;
end

x = -1:0.1:1;
nonl = cdf('norm', x, 0.5,0.2)-0.01;
figure
plot(x, nonl)

% individual subunits
BC = zeros(ncones, sample_size);
for i = 1:ncones
    BC(i, :) = cdf('norm', cones(i,:), 0.5,0.2)-0.01;
%     figure
%     plot(BC(i, :))
end

% summing subunits
BC = zeros(ncones-1, sample_size);
BC(1,:) = cdf('norm', sum(cones(1:2,:)), 0.5,0.2)-0.01;
for i = 3:ncones
    BC(i-1, :) = cdf('norm', cones(i,:), 0.5,0.2)-0.01;
%     figure
%     plot(BC(i, :))
end


% RGC
RGC = sum(BC);
figure
plot(RGC)

hist(RGC)

x1 = -1:0.1:3;
nonl1 = cdf('norm', x1, 1.5,1)-0.06;
figure
plot(x1, nonl1)

res = cdf('norm', RGC, 1.5,1)-0.06;
res(RGC<0.25) = 0;
figure
plot(res)

% calculate conditioned responses
cone1 = 1;
nbins_cone1 = 6;
nbins_cone2 = 6;
size_cone1_bin = floor(size(cones,2)/nbins_cone1);
size_cone2_bin = floor(size(cones,2)/nbins_cone2);
tmp = sort(cones(cone1,:));
cone1_contr = tmp(1:size_cone1_bin:end);

figure
for cone2 = 2:ncones
    cond_rate = zeros(nbins_cone2-1, nbins_cone1-1);

    % second cone: bin the inputs
    tmp = sort(cones(cone2,:));
    cone2_contr = tmp(1:size_cone2_bin:end);
    
    % go through contrast bins on cone2 and calculate the response
    % of cone 1
    cc=1;
    for i=[1,nbins_cone2]
        % find instances when cone2 had certain contrast
        cone2_valid = find(cones(cone2,:)>cone2_contr(i) & cones(cone2,:)<=cone2_contr(i+1));
        
        % go through contrast range on cone1
        for j = 1:nbins_cone1
            cone1_valid = find(cones(cone1,:)>cone1_contr(j) & cones(cone1,:)<=cone1_contr(j+1));
            cone12_valid = intersect(cone1_valid, cone2_valid);
            cond_rate(i,j) = mean(res(cone12_valid));
        end
        subplot(1,2,cc)
        if cone2==2
            plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cond_rate(i,:)-mean(cond_rate(i,:)), 'color', 'b', 'linewidth',2);
        else
            plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cond_rate(i,:)-mean(cond_rate(i,:)));
        end
        hold on
        cc = cc+1;
    end
    
end




