
col = [1 0 0; 0 0 1; 0.7 0 0.7; 0.1 0.5 0.1];
%% NDF5
path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data005/data005';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

figure
hold on
scale = 20;
a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt,tmp, '*', 'color', col(j,:));    
end

%% NDF4

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data008/data008';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);


scale = 10;
a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID, scale);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+2,tmp, '*', 'color', col(j,:));     
end

%% NDF3

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data014/data014';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);


scale = 8;
a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID, scale);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+4,tmp, '*', 'color', col(j,:));     
end



%% NDF2

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data018/data018';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);


scale = 8;
a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID, scale);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+6,tmp, '*', 'color', col(j,:));     
end



%% NDF1

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data022/data022';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);


scale = 8;
a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID, scale);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+8,tmp, '*', 'color', col(j,:));     
end


%% NDF0

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data026/data026';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);


scale = 4;
a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID, scale);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+10,tmp, '*', 'color', col(j,:));     
end

%%
legend('ON parasol', 'OFF parasol', 'ON midget', 'OFF midget')

set(gca, 'xtick', 0.5:2:11, 'xticklabel', {'NDF5', 'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'})

title('2015-03-09-2')



%%






sta = squeeze(datarun.stas.stas{cellID});
tc = datarun.vision.timecourses(cellID).g(11:end);
sta = sta(:,:,11:end);
sta = sta.*reshape(repmat(tc', size(sta,1)^2,1), size(sta,1), size(sta,2),20 );
sta = sum(sta,3);
sta = sta/(max(sta(:))*2)+0.5;
sta = imresize(sta, scale);

figure
colormap gray
imshow(1-sta)
hold on

sta = uint8(round(sta*256));
t = sta([1:30 290:end],[1:30 290:end])*1.2;
tmp = double(max(t(:)))/255;
bw = im2bw(sta, tmp);
%             bw = im2bw(sta, graythresh(sta));
bw2 = imfill(bw,'holes');

%             figure
%             imshow(bw2)

L = bwlabel(bw2);
tmp = [];
for k = 1:length(unique(L(:)))
    tmp(k) = nnz(L==k);
end
[~, myInd] = max(tmp);
L(L~=myInd) = 0;
L(L==myInd) = 1;
[r,c]=find(L==1,1);
contour = bwtraceboundary(L,[r c],'W',8,Inf,'counterclockwise');
hold on;
plot(contour(:,2),contour(:,1),'r','LineWidth',1);
t = nnz(L);

figure
plot_rf_summaries(datarun, datarun.cell_ids(cellID), 'scale', 20,'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')


figure
sta = squeeze(datarun.stas.stas{cellID});
imagesc(imresize(sta(:,:,26), 20, 'method', 'nearest'))
hold on
colormap gray
plot_rf_summaries(datarun, datarun.cell_ids(cellID),'scale', 20,'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')



figure
sta = squeeze(datarun.stas.stas{cellID});
imagesc(sta(:,:,26))
hold on
colormap gray
set(gca, 'dataaspectratio', [1 1 1])
plot_rf_summaries(datarun, datarun.cell_ids(cellID),'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')
axis([0.5 16.5 0.5 16.5])


figure
plot_rf_summaries(datarun, datarun.cell_ids(cellID),'scale', 20,'clear', false, 'label', true, 'label_color', 'k', 'plot_fits', true, 'fit_color', 'k')
hold on
plot(contour(:,2)+scale/2,contour(:,1)+scale/2,'r','LineWidth',1);

t = int8(zeros(320));
t(contour(:,2),contour(:,1)) = 1;
figure
imagesc(t)
a = roipoly(t,contour(:,2)+scale/2,contour(:,1)+scale/2);
figure
imagesc(a)
hold on

