
col = [1 0 0; 0 0 1; 0.7 0 0.7; 0.1 0.5 0.1];
% 512 system (640x320)
% data000 (6.3 BW-32-8), data001 (5.1 BW-16-4)
% data002 (3.9 BW-16-4), data004 (3.9 BW-32-4)
% data005 (3.0 BW-16-4), data007 (2.3 RGB-8-3)
% data009 (2.0 RGB-8-2), data011 (1.6 RGB-8-2)
% data013 (1.3 RGB-8-1), data015 (1.3 BW-16-4)
% data016 (0.6 RGB-8-1), data018 (0.6 BW-3-4 212x106 last 2 bin files deleted)

chrono_order = [0:2:20];
% nd_order = [14 12 10 4 0 16 6 2 8 18];

ord = chrono_order;

%% NDF5.1
path2data = '/Volumes/Analysis/2007-03-02-1/data001-map-gdf/data001-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

figure
hold on
a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(1),tmp, '*', 'color', col(j,:));    
end

%% NDF3.9

path2data = '/Volumes/Analysis/2007-03-02-1/data002-map-gdf/data002-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(2),tmp, '*', 'color', col(j,:));     
end

%% NDF3.9

path2data = '/Volumes/Analysis/2007-03-02-1/data004-map-gdf/data004-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(3),tmp, '*', 'color', col(j,:));     
end



%% NDF3.0

path2data = '/Volumes/Analysis/2007-03-02-1/data005-map-gdf/data005-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(4),tmp, '*', 'color', col(j,:));     
end



%% NDF2.3

path2data = '/Volumes/Analysis/2007-03-02-1/data007-map-gdf/data007-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);


a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(5),tmp, '*', 'color', col(j,:));     
end


%% NDF2.0

path2data = '/Volumes/Analysis/2007-03-02-1/data009-map-gdf/data009-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(6),tmp, '*', 'color', col(j,:));     
end

%% NDF1.6

path2data = '/Volumes/Analysis/2007-03-02-1/data011-map-gdf/data011-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(7),tmp, '*', 'color', col(j,:));     
end

%% NDF1.3

path2data = '/Volumes/Analysis/2007-03-02-1/data013-map-gdf/data013-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(8),tmp, '*', 'color', col(j,:));     
end


%% NDF1.3

path2data = '/Volumes/Analysis/2007-03-02-1/data015-gdf/data015';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(9),tmp, '*', 'color', col(j,:));     
end


%% NDF0.6

path2data = '/Volumes/Analysis/2007-03-02-1/data016-gdf/data016/data016';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(10),tmp, '*', 'color', col(j,:));     
end


%%
legend('ON parasol', 'OFF parasol', 'ON midget', 'OFF midget')

set(gca, 'xtick', 0.5:2:19, 'xticklabel', {'NDF5.1', 'NDF3.9', 'NDF3.9', 'NDF3.0', 'NDF2.3', 'NDF2.0', 'NDF1.6', 'NDF1.3', 'NDF1.3', 'NDF0.6'})
title('2007-03-02-1')




%% NDF3.9
% data002 (3.9 BW-16-4), data004 (3.9 BW-32-4)

path2data = '/Volumes/Analysis/2007-03-02-1/data002-map-gdf/data002-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

figure
for j=1:4
    subplot(2,2,j)
    hold on
    plot_rf_summaries(datarun, {j}, 'scale', 1,'clear', false, 'label', false, 'plot_fits', true, 'fit_color', 'k')     
end

path2data = '/Volumes/Analysis/2007-03-02-1/data004-map-gdf/data004-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

for j=1:4
    subplot(2,2,j)
    hold on
    plot_rf_summaries(datarun, {j}, 'scale', 2,'clear', false, 'label', false, 'plot_fits', true, 'fit_color', 'r')     
end

figure
a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = mean(datarun.vision.sta_fits{cellID}.sd) ;
    end    
    tmp = a(a(:,j)>0,j);
    subplot(2,2,j)
    plot(tmp, 'x-', 'color', col(j,:));   
end

path2data = '/Volumes/Analysis/2007-03-02-1/data002-map-gdf/data002-map';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);
a = zeros(100,4);

for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = mean(datarun.vision.sta_fits{cellID}.sd) ;
    end    
    tmp = a(a(:,j)>0,j);
    subplot(2,2,j)
    hold on
    plot(tmp, '+-', 'color', col(j,:));   
end

title('ON parasol')
title('OFF parasol')
title('ON midget')
title('OFF midget')

