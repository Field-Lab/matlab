%% 2011-10-25-9
%% load stuff

% map
vormap = load('/Volumes/Data/2011-10-25-9/Visual/2011-10-25-9_f06_vorcones/map-0000.txt');
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


vorrun = load_data('/Volumes/Analysis/2011-10-25-9/d02-10-norefit/data008/data008');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_sta(vorrun);
vorrun = load_neurons(vorrun);



vorrun1 = load_data('/Volumes/Analysis/2011-10-25-9/d02-10-norefit/data009/data009');
vorrun1 = load_params(vorrun1,'verbose',1);
vorrun1 = load_sta(vorrun1);
vorrun1 = load_neurons(vorrun1);


cnt = 1

for i=vorrun.cell_types{1}.cell_ids  
    cellID = find(vorrun.cell_ids==i);
    
    if cnt==26
        cnt=1;figure;
    end

    sta = squeeze(vorrun.stas.stas{cellID});
    if ~any(isnan(sta(:)))
        subplot(5,5,cnt)
        cnt = cnt+1;
%         sta = -sta;
        [mm, frame] = find(sta==max(sta(:)),1);
        a = find(sta(:,frame)>sta(mm,frame)*0.4);
        mean_sta = mean(sta(a,:));
        sta = sta.*repmat(mean_sta,1264,1);
        sta = sum(sta,2);
        sta = sta/max(sta);
        [t,ic] = sort(sta);
        plot(t)
        hold on
        
        sta = squeeze(vorrun1.stas.stas{cellID});
%         sta = -sta;
%         [mm, frame] = find(sta==max(sta(:)),1);
        a = find(sta(:,frame)>sta(mm,frame)*0.4);
        mean_sta = mean(sta(a,:));
        sta = sta.*repmat(mean_sta,1264,1);
        sta = sum(sta,2);
        sta = sta/max(sta);
        [t,ic] = sort(sta);
        plot(t)
        axis([1200,1270,0,Inf])
        title(['datarun ID ', int2str(cellID), ', Cell ID ', int2str(i)])
        drawnow
    end
end
    

% sta direct comparison
for cellID=1: length(vorrun.cell_ids ) 
    datarunID = vorrun.cell_ids(cellID);
    
    sta = squeeze(vorrun.stas.stas{cellID});
    sta1 = squeeze(vorrun1.stas.stas{cellID});
    if ~any(isnan(sta(:))) & ~any(isnan(sta1(:)))
        
        figure
        set(gcf, 'position', [-1714         532        1300         561])
        [mm, frame] = find(abs(sta)==max(abs(sta(:))),1);
        if max(abs(sta(:)))==max(sta(:))
            a = find(sta(:,frame)>sta(mm,frame)*0.4);
            pol = 1;
        else
            a = find(sta(:,frame)<sta(mm,frame)*0.4);
            pol = -1;
        end
        if length(a)==1
            mean_sta = sta(a,:);
        else
            mean_sta = mean(sta(a,:));
        end
        sta = sta.*repmat(pol*mean_sta,1264,1);
        sta = sum(sta,2);
        [full_sta, ~] = expand_voronoi_sta(sta, vormap);
        subplot(1,2,1)
        colormap gray
        imagesc(full_sta)
        set(gca, 'dataaspectratio', [1 1 1])
        title(['datarun ID ', int2str(cellID), ', Cell ID ', int2str(datarunID), ', BW-3 NORMAL'])
        

        sta = sta1;
        [mm, frame] = find(abs(sta)==max(abs(sta(:))),1);
        if max(abs(sta(:)))==max(sta(:))
            a = find(sta(:,frame)>sta(mm,frame)*0.4);
            pol = 1;
        else
            a = find(sta(:,frame)<sta(mm,frame)*0.4);
            pol = -1;
        end
        if length(a)==1
            mean_sta = sta(a,:);
        else
            mean_sta = mean(sta(a,:));
        end
        sta = sta.*repmat(pol*mean_sta,1264,1);
        sta = sum(sta,2);
        [full_sta, ~] = expand_voronoi_sta(sta, vormap);
        subplot(1,2,2)
        colormap gray
        imagesc(full_sta)
        set(gca, 'dataaspectratio', [1 1 1])

        title(['datarun ID ', int2str(cellID), ', Cell ID ', int2str(datarunID), ', BW-3, SPARSE 0.5'])
        drawnow
        
        saveas(gcf, ['/Users/alexth/Desktop/SPARSE_Voronoi/2011-10-25-9/svg/Cell_', int2str(cellID), '_datarunID_', int2str(datarunID), '.svg'])
        saveas(gcf, ['/Users/alexth/Desktop/SPARSE_Voronoi/2011-10-25-9/tiff/Cell_', int2str(cellID), '_datarunID_', int2str(datarunID), '.tiff'])
        
        close all
    end
end
    
%% 2011-12-13-2
%% load stuff
clear
% map
vormap = load('/Volumes/Data/2011-12-13-2/Visual/2011-12-13-2_f14_vorcones/map-0000.txt');
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


vorrun = load_data('/Volumes/Analysis/2011-12-13-2/d04-20-norefit/data017-from-d04_20/data017-from-d04_20');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_sta(vorrun);
vorrun = load_neurons(vorrun);



vorrun1 = load_data('/Volumes/Analysis/2011-12-13-2/d04-20-norefit/data019-from-d04_20/data019-from-d04_20');
vorrun1 = load_params(vorrun1,'verbose',1);
vorrun1 = load_sta(vorrun1);
vorrun1 = load_neurons(vorrun1);

 

cnt = 1

for i=vorrun.cell_types{3}.cell_ids  
    cellID = find(vorrun.cell_ids==i);
    
    if cnt==26
        cnt=1;figure;
    end

    sta = squeeze(vorrun.stas.stas{cellID});
    if ~any(isnan(sta(:)))
        subplot(5,5,cnt)
        cnt = cnt+1;
%         sta = -sta;
[mm, frame] = find(sta==max(sta(:)),1);
a = find(sta(:,frame)>sta(mm,frame)*0.4);
if length(a)==1
    mean_sta = sta(a,:);
else
    mean_sta = mean(sta(a,:));
end
        sta = sta.*repmat(mean_sta,954,1);
        sta = sum(sta,2);
        sta = sta/max(sta);
        [t,ic] = sort(sta);
        plot(t)
        hold on
        
        sta = squeeze(vorrun1.stas.stas{cellID});
%         sta = -sta;
        [mm, frame] = find(sta==max(sta(:)),1);
        a = find(sta(:,frame)>sta(mm,frame)*0.4);
        if length(a)==1
            mean_sta = sta(a,:);
        else
            mean_sta = mean(sta(a,:));
        end
        sta = sta.*repmat(mean_sta,954,1);
        sta = sum(sta,2);
        sta = sta/max(sta);
        [t,ic] = sort(sta);
        plot(t)
        axis([860,960,0,Inf])
        title(['datarun ID ', int2str(cellID), ', Cell ID ', int2str(i)])
        drawnow
    end
end
    

% sta direct comparison
for cellID=1: length(vorrun.cell_ids ) 
    datarunID = vorrun.cell_ids(cellID);
    
    sta = squeeze(vorrun1.stas.stas{cellID});
    sta1 = squeeze(vorrun.stas.stas{cellID});
    if ~any(isnan(sta(:))) & ~any(isnan(sta1(:)))
        
        figure
        set(gcf, 'position', [-1714         532        1300         561])
        [mm, frame] = find(abs(sta)==max(abs(sta(:))),1);
        if max(abs(sta(:)))==max(sta(:))
            a = find(sta(:,frame)>sta(mm,frame)*0.4);
            pol = 1;
        else
            a = find(sta(:,frame)<sta(mm,frame)*0.4);
            pol = -1;
        end
        if length(a)==1
            mean_sta = sta(a,:);
        else
            mean_sta = mean(sta(a,:));
        end
        sta = sta.*repmat(pol*mean_sta,954,1);
        sta = sum(sta,2);
        [full_sta, ~] = expand_voronoi_sta(sta, vormap);
        subplot(1,2,1)
        colormap gray
        imagesc(full_sta)
        set(gca, 'dataaspectratio', [1 1 1])
        title(['datarun ID ', int2str(cellID), ', Cell ID ', int2str(datarunID), ', BW-6 NORMAL'])
        

        sta = sta1;
        [mm, frame] = find(abs(sta)==max(abs(sta(:))),1);
        if max(abs(sta(:)))==max(sta(:))
            a = find(sta(:,frame)>sta(mm,frame)*0.4);
            pol = 1;
        else
            a = find(sta(:,frame)<sta(mm,frame)*0.4);
            pol = -1;
        end
        if length(a)==1
            mean_sta = sta(a,:);
        else
            mean_sta = mean(sta(a,:));
        end
        sta = sta.*repmat(pol*mean_sta,954,1);
        sta = sum(sta,2);
        [full_sta, ~] = expand_voronoi_sta(sta, vormap);
        subplot(1,2,2)
        colormap gray
        imagesc(full_sta)
        set(gca, 'dataaspectratio', [1 1 1])

        title(['datarun ID ', int2str(cellID), ', Cell ID ', int2str(datarunID), ', BW-6, SPARSE 0.5'])
        drawnow
        
        saveas(gcf, ['/Users/alexth/Desktop/SPARSE_Voronoi/2011-12-13-2/svg/Cell_', int2str(cellID), '_datarunID_', int2str(datarunID), '.svg'])
        saveas(gcf, ['/Users/alexth/Desktop/SPARSE_Voronoi/2011-12-13-2/tiff/Cell_', int2str(cellID), '_datarunID_', int2str(datarunID), '.tiff'])
        
        close all
    end
end
