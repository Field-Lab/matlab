%% load stuff

local_path = '/Volumes/Analysis/';

vormap = load([local_path, '2016-02-17-4/stimuli/maps/map_data001_from_data000_600_1500.txt']);

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
        [r, c] = find(dd,1);
        contour = bwtraceboundary(dd,[r c],'W',8,Inf,'counterclockwise');
        contour= round(contour/5)+0.5;
        plot(contour(:,2), contour(:,1), 'linewidth', 2)
        voronoi_contours{i, 1} = contour(:,2);
        voronoi_contours{i, 2} = contour(:,1);
    end
end

%Coarse run run 1
datarunc = load_data([local_path, '2016-02-17-4/d00-05-norefit/data000/data000']);
datarunc = load_params(datarunc,'verbose',1);
datarunc = load_sta(datarunc);


%SC run 1
datarun = load_data([local_path, '2016-02-17-4/d00-05-norefit/data001/data001']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);


%SC run OTHER DATE
datarun7 = load_data([local_path, '2016-02-17-7/d01-13-norefit/data002/data002']);
datarun7 = load_params(datarun7,'verbose',1);
datarun7 = load_sta(datarun7);


%SC run 1 FALSE
datarunf = load_data([local_path, '2016-02-17-4/d00-05-norefit/data001_false/data001']);
datarunf = load_params(datarunf,'verbose',1);
datarunf = load_sta(datarunf);

%SC run 2
datarun2 = load_data([local_path, '2016-02-17-4/d00-05-norefit/data005/data005']);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2);

% voronoi run
vorrun = load_data([local_path, '2016-02-17-4/d00-05-norefit/data004/data004']);
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_sta(vorrun);
vorrun = load_neurons(vorrun);
[inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-1-0.48-11111-1340x1-60.35.xml');
inputs = squeeze(inputs);


a = load('/Volumes/Analysis/2016-02-17-4/cone_data/manual/map_data001_manual_postexp_info.mat')
b = load('/Volumes/Analysis/2016-02-17-4/cone_data/manual/map_data005_manual_postexp_info.mat')
map = load('/Volumes/Analysis/2016-02-17-4/stimuli/maps/map_data001_from_data000_600_1500.txt');
figure
imagesc(map)
hold on
plot(b.cones(:,1), b.cones(:,2), 'r+');
plot(a.cones(:,1), a.cones(:,2), 'bx');
cones1 = a.cones;
cones2 = b.cones; 

%% cone profiles: after
good_cells = [541, 556, 571, 647, 796, 1276, 1367, 1486, 1816, 1892, 1921, 2087, 2117,...
    2161, 2192, 2222, 2386, 2658, 3001, 3031, 3046, 3182, 3258, 3272, 3721, 4111, 4486,...
    4966, 4967, 5058, 5341, 5613, 5822, 5914, ...
    6017, 6061, 6244, 6256, 6332, 6483, 6572, 6991, 7156, 7246, 7502, 7652, 7667, 7697];
all_sta = zeros(600,800,length(good_cells));
cnt = 1;
pols = zeros(1,length(good_cells));
max_pics = pols;
for i = 1:length(datarun.cell_ids)
    if ~isempty(find(good_cells == datarun.cell_ids(i),1))
        tmp = datarun.stas.stas{i};
        t = zeros(1,6);
        for j=1:6
            sta=squeeze(tmp(:,:,:,j));
            t(j) = max(abs(sta(:)));
        end
        [~, frame] = max(t);
        sta=tmp(:,:,:,frame);
        if max(sta(:))<max(abs(sta(:))) % OFF cell
            pols(cnt) = -1;
        end
        all_sta(:,:,cnt)=imresize(double(squeeze(sta)), 2, 'method', 'nearest');
        tmp = sta*pols(cnt);
        max_pics(cnt) = max(tmp(:));
        cnt = cnt+1;
    end
end

cone_portrait = cell(max(vormap(:)),2);
for i=1:max(vormap(:))
    
    [a,b] = find(vormap==i);
    tmp = 0;
    for j=1:length(a)
        tmp = tmp+squeeze(all_sta(a(j),b(j),:));
    end
    tmp = tmp'.*pols;
    cone_portrait{i,1} = find(tmp>(max(tmp)*0.85) & tmp>max_pics*0.9);
    cone_portrait{i,2} = [a b];
end


bord = 2;
cnt = 1;
combo = zeros(600,800);
for i=1:max(vormap(:))
    tmp = cone_portrait{i,1};
    if ~isempty(tmp)
        tt = 0;
        for j=1:length(tmp)
            tt = tt+pols(tmp(j))*all_sta(:,:,tmp(j));
        end
        tmp =  cone_portrait{i,2};
        combo(min(tmp(:,1))-bord:max(tmp(:,1))+bord, min(tmp(:,2))-bord:max(tmp(:,2))+bord) = tt(min(tmp(:,1))-bord:max(tmp(:,1))+bord, min(tmp(:,2))-bord:max(tmp(:,2))+bord);
    end
    cnt = cnt+1;
end
figure
imagesc(combo)
hold on
for i=1:max(vormap(:))
    plot(voronoi_contours{i, 1}, voronoi_contours{i, 2}, 'color', 'r')
end


for k = 1:9:100
    
    figure
    set(gcf, 'position', [-1580         125        1085         966])
    bord = 15;
    cnt = 1;
    for i=k:k+8
        tmp = cone_portrait{i,1};
        if ~isempty(tmp)
            tt = 0;
            for j=1:length(tmp)
                tt = tt+pols(tmp(j))*all_sta(:,:,tmp(j));
            end
            subplot(3,3,cnt)
            imagesc(tt)
            tmp =  cone_portrait{i,2};
            axis([min(tmp(:,2))-bord  max(tmp(:,2))+bord min(tmp(:,1))-bord  max(tmp(:,1))+bord ])
            hold on
            plot(voronoi_contours{i, 1}, voronoi_contours{i, 2}, 'color', 'r')
            title(int2str(j))
        end
        cnt = cnt+1;
    end
    
end

  


%% plot stuff
save_path = [local_path, '2016-02-17-4/data004_d00-05/'];
nbins_cone1 = 8;
nbins_cone2 = 8;
sta_params.length = 20;
sta_params.offset = 0;
fraction = 0.9;

% denoise 
tmp = 0;
for i=1:length(datarun.cell_ids)
    tmp = tmp + squeeze(datarun.stas.stas{i}(:,:,:,4));
end
noise_sta = tmp/length(datarun.cell_ids);
tmp = 0;
for i=1:length(datarun2.cell_ids)
    tmp = tmp + squeeze(datarun2.stas.stas{i}(:,:,:,4));
end
noise_sta2 = tmp/length(datarun2.cell_ids);

tmp = 0;
for i=1:length(vorrun.cell_ids)
    tmp = tmp + squeeze(vorrun.stas.stas{i});
end
noise_vor = repmat(mean(tmp/length(vorrun.cell_ids),2), 1,30);

        
[~, cones] = expand_voronoi_sta(noise_vor, vormap);
        

for kkk= vorrun.cell_types{4}.cell_ids(1:end)  %vorrun.cell_ids(22:end)
    close all
    datarunID = find(vorrun.cell_ids==kkk);
    
    for i=1:length(vorrun.cell_types)
        if ~isempty(find(vorrun.cell_types{i}.cell_ids==kkk, 1))
            cell_type = vorrun.cell_types{i}.name;
            break
        end
    end
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID})-noise_vor;
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
        [~,pos] = find(raw_sta==b,1);
        a = b;
    end
    tmp = raw_sta * pol;
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    
    if abs(a) > thresh*1.5 && length(spikes)>300

        center_cones = find(sum(tmp(:,pos-1:pos+1)>thresh,2));
        center_cones(isnan(cones(center_cones,1)))=[];
        
        if length(center_cones)>1
            
            
            %     figure
            %     plot(raw_sta(center_cones,:)')
            timecourse = mean(tmp(center_cones,10:end));
%                 figure
%                 plot(timecourse)
            
%             if length(center_cones)>30
%                 [~, a] = sort(conv_sta);
%                 center_cones = sort(a(end-29:end));
%             end
            conv_sta = sum(tmp(:,10:end).*repmat(timecourse,size(tmp,1),1),2);
            
            % prepare voronoi and single cone STAs
            [full_sta, ~] = expand_voronoi_sta(conv_sta, vormap);
            full_sta = pol*full_sta/(max(full_sta(:))*1.5);
            sta1 = squeeze(datarun.stas.stas{datarunID}(:,:,:,4))-noise_sta;
            sta1 = imresize(sta1, 2, 'method', 'nearest');
            sta2 = squeeze(datarun2.stas.stas{datarunID}(:,:,:,4))-noise_sta2;
            sta2 = imresize(sta2, 2, 'method', 'nearest');
            
            raw_sta = raw_sta(center_cones, end-1:-1:10);
            
            
            filt_inputs = zeros(length(center_cones), size(inputs,2)-sta_params.length+1);
            cnt = 1;
            for current_cone=center_cones'
                filt_inputs(cnt,:)=conv(inputs(current_cone,:), raw_sta(cnt,:),'valid');
                cnt=cnt+1;
            end
            spikes_tmp = spikes;
            spikes_tmp(spikes<sta_params.length) = [];
            
            offline_contrast_response(filt_inputs, spikes_tmp-sta_params.length+1,...
                center_cones, vormap, cones, full_sta, datarunID,...
                kkk, [save_path int2str(nbins_cone2) '_bins/', cell_type '/'], sta1, sta2,...
                cones1, cones2);
            
            tic
            [resnorm_all, resnorm_all1] = cone_interaction_fits(filt_inputs, spikes_tmp-sta_params.length+1,...
                nbins_cone1, nbins_cone2, center_cones, vormap, cones, full_sta, datarunID,...
                kkk, [save_path int2str(nbins_cone2) '_bins/', cell_type '/'], sta1, sta2,...
                cones1, cones2, voronoi_contours, pol);
            toc
            
            [a, b] = cone_interaction_fits(filt_inputs, spikes_tmp-sta_params.length+1,...
                nbins_cone1, nbins_cone2, center_cones, vormap, cones, full_sta, datarunID,...
                kkk, [save_path int2str(nbins_cone2) '_bins/', cell_type '/'], sta1, sta2,...
                cones1, cones2, voronoi_contours, pol);
            
            figure
            for coneid=1:20
                subplot(4,5,coneid)
                plot(resnorm_all(coneid,:), resnorm_all1(coneid,:), 'o', 'markersize', 15)
                hold on
                tmp = max([resnorm_all(coneid,:) resnorm_all1(coneid,:)]);
                line([0,tmp], [0 tmp], 'color', 'k')
                xlabel('x shift')
                ylabel('y shift')
                axis([0 tmp 0 tmp])
%                 axis([0 0.1 0 0.1])
            end
            
            figure
            for coneid=21:40
                subplot(4,5,coneid-20)
                plot(resnorm_all(coneid,:), resnorm_all1(coneid,:), 'o', 'markersize', 15)
                hold on
                tmp = max([resnorm_all(coneid,:) resnorm_all1(coneid,:)]);
                line([0,tmp], [0 tmp], 'color', 'k')
                xlabel('x shift')
                ylabel('y shift')
                axis([0 tmp 0 tmp])
%                 axis([0 0.1 0 0.1])
            end
            
            figure
            for coneid=41:60
                subplot(4,5,coneid-40)
                plot(resnorm_all(coneid,:), resnorm_all1(coneid,:), 'o', 'markersize', 15)
                hold on
                tmp = max([resnorm_all(coneid,:) resnorm_all1(coneid,:)]);
                line([0,tmp], [0 tmp], 'color', 'k')
                xlabel('x shift')
                ylabel('y shift')
                axis([0 tmp 0 tmp])
            end
            
            figure
            for coneid=61:80
                subplot(4,5,coneid-60)
                plot(resnorm_all(coneid,:), resnorm_all1(coneid,:), 'o', 'markersize', 15)
                hold on
                tmp = max([resnorm_all(coneid,:) resnorm_all1(coneid,:)]);
                line([0,tmp], [0 tmp], 'color', 'k')
                xlabel('x shift')
                ylabel('y shift')
                axis([0 tmp 0 tmp])
            end
                
        end
    end
end


figure
hold on
axis ij
for i=1:length(center_cones)
    plot(voronoi_contours{center_cones(i),1}, voronoi_contours{center_cones(i),2})
end

coneid =5;
figure
plot(resnorm_all(coneid,:), resnorm_all1(coneid,:), 'o', 'markersize', 15)
              
t = find(resnorm_all1(coneid,:)>0.01);
t1 = find(resnorm_all1(coneid,:)<0.01);
figure
hold on
axis ij
for i=t
    plot(voronoi_contours{center_cones(i),1}, voronoi_contours{center_cones(i),2}, 'r')
end
for i=t1
    plot(voronoi_contours{center_cones(i),1}, voronoi_contours{center_cones(i),2}, 'b')
end
plot(voronoi_contours{center_cones(coneid),1}, voronoi_contours{center_cones(coneid),2}, 'k', 'linewidth',3)


for coneid=1:20

           
t = find(resnorm_all1(coneid,:)>0.006);
t1 = find(resnorm_all1(coneid,:)<0.006);
tt(coneid) = length(t);
tt1(coneid) = length(t1);
end

coneid =18;
figure
hold on
axis ij
for i=t
    plot(voronoi_contours{center_cones(i),1}, voronoi_contours{center_cones(i),2}, 'r')
end
for i=t1
    plot(voronoi_contours{center_cones(i),1}, voronoi_contours{center_cones(i),2}, 'b')
end
plot(voronoi_contours{center_cones(coneid),1}, voronoi_contours{center_cones(coneid),2}, 'k', 'linewidth',3)



figure
for coneid=1:20
    subplot(4,5,coneid)
    tmp = resnorm_all{coneid}(2:end);
    tmp1 = resnorm_all1{coneid}(2:end);
    plot(tmp(1:2:end), tmp1(1:2:end), 'ob', 'markersize', 15)
    hold on
    plot(tmp(2:2:end), tmp1(2:2:end), 'or', 'markersize', 15)
    plot(resnorm_all{coneid}(1), resnorm_all1{coneid}(1), 'xk', 'markersize', 15)
    plot(resnorm_all{coneid}(1), resnorm_all1{coneid}(1), 'ok', 'markersize', 15)
    tmp = max([tmp; tmp1]);
    line([0,tmp], [0 tmp], 'color', 'k')
    xlabel('x shift')
    ylabel('y shift')
    axis([0 tmp 0 tmp])
    title(int2str(coneid))
%     axis([0 0.1 0 0.1])
end


figure
for coneid=21:40
    subplot(4,5,coneid-20)
    tmp = resnorm_all{coneid}(2:end);
    tmp1 = resnorm_all1{coneid}(2:end);
    plot(tmp(1:2:end), tmp1(1:2:end), 'ob', 'markersize', 15)
    hold on
    plot(tmp(2:2:end), tmp1(2:2:end), 'or', 'markersize', 15)
    plot(resnorm_all{coneid}(1), resnorm_all1{coneid}(1), 'xk', 'markersize', 15)
    plot(resnorm_all{coneid}(1), resnorm_all1{coneid}(1), 'ok', 'markersize', 15)
    tmp = max([tmp; tmp1]);
    line([0,tmp], [0 tmp], 'color', 'k')
    xlabel('x shift')
    ylabel('y shift')
    axis([0 tmp 0 tmp])
    title(int2str(coneid))
%     axis([0 0.1 0 0.1])
end



