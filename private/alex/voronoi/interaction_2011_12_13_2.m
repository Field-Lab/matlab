%% load stuff

local_path = '/Volumes/Analysis/';

vormap = load('/Volumes/Data/2011-12-13-2/Visual/2011-12-13-2_f04_vorcones/map-0000.txt');

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

%SC run 1
datarun = load_data([local_path, '2011-12-13-2/d08-11-norefit/data008/data008']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);

%SC run 2
datarun2 = load_data([local_path, '2011-12-13-2/d08-11-norefit/data011/data011']);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2);

% voronoi run
vorrun = load_data([local_path, '2011-12-13-2/d08-11-norefit/data009/data009']);
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_sta(vorrun);
vorrun = load_neurons(vorrun);
[inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-937x1-60.35.xml');
inputs = squeeze(inputs);


a = load('/Volumes/Analysis/2011-12-13-2/cone_data/manual/map_data008_manual_postexp_info.mat')
b = load('/Volumes/Analysis/2011-12-13-2/cone_data/manual/map_data011_manual_postexp_info.mat')
map = load('/Volumes/Data/2011-12-13-2/Visual/2011-12-13-2_f04_vorcones/map-0000.txt');
figure
imagesc(map)
hold on
plot(b.cones(:,1), b.cones(:,2), 'r+');
plot(a.cones(:,1), a.cones(:,2), 'bx');
cones1 = a.cones;
cones2 = b.cones; 


%% plot stuff
save_path = [local_path, '2011-12-13-2/data009_d08-11/'];

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


cone_table = zeros(length(vorrun.cell_ids), length(cones));

% cone lists per cell
cnt = 1;
for kkk= vorrun.cell_ids
    datarunID = find(vorrun.cell_ids==kkk);
    
    for i=1:length(vorrun.cell_types)
        if ~isempty(find(vorrun.cell_types{i}.cell_ids==kkk, 1))
            all_cell_type(cnt) = i;
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
    
    if abs(a) > thresh*1.5
        
        center_cones = find(sum(tmp(:,pos-1:pos+1)>max(tmp(:))*0.3,2));
        
        center_cones(isnan(cones(center_cones,1)))=[];
        
        cone_table(cnt, center_cones) = 1;
    end
    cnt = cnt+1;
end

cone_cells = find(sum(cone_table(:,[549 550]),2)==2);
all_cell_type(cone_cells);


all_cells = [];
for i=1:5
    all_cells = [all_cells vorrun.cell_types{i}.cell_ids(1:end)];
end

for kkk=  [6577 6288] %all_cells(1:end)

    kkk
%     close all
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
        
%         center_cones = find(sum(tmp(:,pos-1:pos+1)>thresh,2));
        center_cones = find(sum(tmp(:,pos-1:pos+1)>max(tmp(:))*0.5,2));
        
        center_cones(isnan(cones(center_cones,1)))=[];
        
        
        if length(center_cones)>1
            length(center_cones)
%             path2save = ['/Volumes/Analysis/2011-12-13-2/data009_d08-11/5_bins/', cell_type, '/'];
            
            
%             timecourse = mean(tmp(center_cones,10:end));
            %                 figure
            %                 plot(timecourse)
            
            
%             conv_sta = sum(tmp(:,10:end).*repmat(timecourse,size(tmp,1),1),2);
%             
%             % prepare voronoi and single cone STAs
%             [full_sta, ~] = expand_voronoi_sta(conv_sta, vormap);
%             full_sta = pol*full_sta/(max(full_sta(:))*1.5);
%             sta1 = squeeze(datarun.stas.stas{datarunID}(:,:,:,4))-noise_sta;
%             sta1 = imresize(sta1, 2, 'method', 'nearest');
%             sta2 = squeeze(datarun2.stas.stas{datarunID}(:,:,:,4))-noise_sta2;
%             sta2 = imresize(sta2, 2, 'method', 'nearest');
% plot_maps_for_interactions(center_cones, vormap, cones, full_sta, datarunID, save_path, sta1, sta2, cones1, cones2, voronoi_contours)
%             
            raw_sta = raw_sta(center_cones, end-1:-1:10);
            
            
            filt_inputs = zeros(length(center_cones), size(inputs,2)-sta_params.length+1);
            cnt = 1;
            for current_cone=center_cones'
                filt_inputs(cnt,:)=conv(inputs(current_cone,:)-noise_vor(current_cone,1), raw_sta(cnt,:),'valid');
                cnt=cnt+1;
            end
            spikes_tmp = spikes;
            spikes_tmp(spikes<sta_params.length) = [];
            
%             err =[];
%             try  [params_ref, params_x, params_y, params_xy,resnorm_x, resnorm_y,resnorm_xy, x, y] = ...
%                 cone_interaction_fits(filt_inputs, spikes_tmp-sta_params.length+1, center_cones);
%             catch err
%             end
            
            
            err =[];
            try  loglikratio = fit_normal_cdfs(filt_inputs, spikes_tmp-sta_params.length+1, center_cones);
            catch err
            end            
            
            if isempty(err)
                save(['/Volumes/Analysis/2011-12-13-2/cone_data/manual/cell_', int2str(kkk), '.mat'], 'loglikratio');
            end
%             
%             if isempty(err)
%                 plot_cone_interaction_fits(params_ref, params_x, params_y, ...
%                     params_xy, resnorm_x, resnorm_y, resnorm_xy, x, y, ...
%                     center_cones, path2save, datarunID, kkk)
%                 
%                 
%                 params_ref_all{cnt1} = params_ref;
%                 params_x_all{cnt1} = params_x;
%                 params_y_all{cnt1} = params_y;
%                 params_xy_all{cnt1} = params_xy;
%                 resnorm_x_all{cnt1} = resnorm_x;
%                 resnorm_y_all{cnt1} = resnorm_y;
%                 resnorm_xy_all{cnt1} = resnorm_xy;
%                 x_all{cnt1} = x;
%                 y_all{cnt1} = y;
%                 center_cones_all{cnt1} = center_cones;
%             end

           
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

% 
% 
% tmp = loglikratio(7, 2, :);
% tmp=reshape(tmp,3,3);
% mean(tmp(:))/std(tmp(:))



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












k = find(all_cells==7351);
center_cones = center_cones_acc{k};
resnorm_all = acc{k};
resnorm_all1 = acc1{k};
save('/Volumes/Analysis/2011-12-13-2/cone_data/manual/fit_data', 'params_ref_all', ...
    'params_x_all', 'params_y_all', 'params_xy_all', 'resnorm_x_all', ...
    'resnorm_y_all', 'resnorm_xy_all', 'x', 'y', 'center_cones_acc', 'all_cells')


for cellid = 24:length( all_cells)
    myCell = all_cells(cellid);
    if ~isempty(acc{cellid})
        datarun_ID = find(vorrun.cell_ids==myCell);
        for i=1:5
            if ~isempty(find(vorrun.cell_types{i}.cell_ids==myCell, 1))
                break
            end
        end
        cell_type = vorrun.cell_types{i}.name;
        resnorm_all = acc{cellid};
        resnorm_all1 = acc1{cellid};
        center_cones = center_cones_acc{cellid};
        clear r c
        for i=1:ceil(length(center_cones)/20)
            if length(center_cones)>20*i
                k = 20;
            else
                k = mod(length(center_cones), 20);
            end
            [r(i), c(i)] = opt_subplots(k);
        end
        path2save = ['/Volumes/Analysis/2011-12-13-2/data009_d08-11/8_bins/', cell_type, '/'];
        coneid = 1;
        for i=1:length(r)
            if length(center_cones)>20*i
                k = 20;
            else
                k = mod(length(center_cones), 20);
            end
            
            figure
            set(gcf, 'position', [-1732         167        1151         938])
            cnt = 1;
            for tmp_cone = coneid:coneid+k-1
                subplot(r(i),c(i),cnt)
                plot(resnorm_all(tmp_cone,:)*100, resnorm_all1(tmp_cone,:)*100, 'o', 'markersize', 5)
                hold on
                tmp = max([resnorm_all(tmp_cone,:)*100 resnorm_all1(tmp_cone,:)*100]);
                line([0,tmp], [0 tmp], 'color', 'k')
                set(gca, 'dataaspectratio', [1 1 1])
                xlabel('x shift', 'Fontsize', 8)
                ylabel('y shift', 'Fontsize', 8)
                axis([0 tmp 0 tmp])
                title([int2str(tmp_cone), ', cone id ', int2str(center_cones(tmp_cone))], 'Fontsize', 10)
                cnt=cnt+1;
            end
            saveas(gcf,[path2save,int2str(datarun_ID), '/tiff/ID_',int2str(myCell),'_cones_', int2str(coneid),'.tiff'])
            saveas(gcf,[path2save,int2str(datarun_ID), '/svg/ID_',int2str(myCell),'_cones_', int2str(coneid),'.svg'])
            drawnow
            close(gcf)
            coneid = coneid+k;
        end
    end
end

coneid = 1;
figure
set(gcf, 'position', [-1732         167        1151         938])
for tmp_cone = 1:13
    subplot(4,4,tmp_cone)
    plot(resnorm_all(tmp_cone,:), resnorm_all1(tmp_cone,:), 'o', 'markersize', 5)
    hold on
    tmp = max([resnorm_all(tmp_cone,:) resnorm_all1(tmp_cone,:)]);
    line([0,tmp], [0 tmp], 'color', 'k')
    set(gca, 'dataaspectratio', [1 1 1])
    xlabel('x shift')
    ylabel('y shift')
    axis([0 tmp 0 tmp])
    title([int2str(tmp_cone), ', cone id ', int2str(center_cones(tmp_cone))])
end
saveas(gcf,[path2save,int2str(datarun_ID), '/tiff/ID_',int2str(myCell),'_cones_', int2str(coneid),'.tiff'])
saveas(gcf,[path2save,int2str(datarun_ID), '/svg/ID_',int2str(myCell),'_cones_', int2str(coneid),'.svg'])
drawnow
close(gcf)



figure
hold on
axis ij
for i=1:length(center_cones)
    plot(voronoi_contours{center_cones(i),1}, voronoi_contours{center_cones(i),2})
end

coneid =10;
figure
plot(resnorm_all(coneid,:), resnorm_all1(coneid,:), 'o', 'markersize', 15)
              
t = find(resnorm_all1(coneid,:)>0.14);
t1 = find(resnorm_all1(coneid,:)<0.14);
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

