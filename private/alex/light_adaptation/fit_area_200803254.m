
col = [1 0 0; 0 0 1; 0.7 0 0.7; 0.1 0.5 0.1];
% 512 system (640x320)
% NDFs: data000 (0, RGB-10-2), data001 (0, BW-20-4), 
% data003 (0.0 RGB-1-6)

% data004 (0 RGB-8-2), data005 (0.6, RGB-8-2)
% data006 (1.0 RGB-8-2), data007 (1.6 RGB-8-4)
% data008 (3.3 RGB-16-6), data009 (0 RGB-8-2)
% data010 (1.3 RGB-8-3), data011 (2.6 RGB-16-4),
% data012 (1.3 RGB-8-3), data013 (0 RGB-8-4)

chrono_order = [0:2:18];
nd_order = [14 12 10 4 0 16 6 2 8 18];

ord = nd_order;

%% NDF0
path2data = '/Volumes/Analysis/2008-03-25-4/d03-13-norefit/data004/data004';
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

%% NDF0.6

path2data = '/Volumes/Analysis/2008-03-25-4/d03-13-norefit/data005/data005';
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

%% NDF1.0

path2data = '/Volumes/Analysis/2008-03-25-4/d03-13-norefit/data006/data006';
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



%% NDF1.6

path2data = '/Volumes/Analysis/2008-03-25-4/d03-13-norefit/data007/data007';
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



%% NDF3.3

path2data = '/Volumes/Analysis/2008-03-25-4/d03-13-norefit/data008/data008';
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


%% NDF0

path2data = '/Volumes/Analysis/2008-03-25-4/d03-13-norefit/data009/data009';
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

%% NDF1.3

path2data = '/Volumes/Analysis/2008-03-25-4/d03-13-norefit/data010/data010';
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

%% NDF2.6

path2data = '/Volumes/Analysis/2008-03-25-4/d03-13-norefit/data011/data011';
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

path2data = '/Volumes/Analysis/2008-03-25-4/d03-13-norefit/data012/data012';
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


%% NDF0

path2data = '/Volumes/Analysis/2008-03-25-4/d03-13-norefit/data013/data013';
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

% chrono
set(gca, 'xtick', 0.5:2:19, 'xticklabel', {'NDF0', 'NDF0.6', 'NDF1.0', 'NDF1.6', 'NDF3.3', 'NDF0', 'NDF1.3', 'NDF2.6', 'NDF1.3', 'NDF0'})
title('2007-02-06-4, chrono order')

% ndf
set(gca, 'xtick', 0.5:2:19, 'xticklabel', {'NDF3.3', 'NDF2.6', 'NDF1.6', 'NDF1.3', 'NDF1.3', 'NDF1.0', 'NDF0.6', 'NDF0', 'NDF0', 'NDF0'})
title('2007-02-06-4, NDF order')

