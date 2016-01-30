
col = [1 0 0; 0 0 1; 0.7 0 0.7; 0.1 0.5 0.1];
% 512 system (640x480)
% NDFs: data001 (5.1, BW-16-4), data004 (3.9 BW-16-4), 
% data005 (3.0 BW-16-4), data007 (2.3 RGB-16-3), data009 (2, RGB-8-2)
% data011 (1.6 RGB-8-2), data013 (1.3 RGB-8-1), data015 (0.6 RGB-8-1)

%% NDF5.1
path2data = '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/data001/data001';
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
    plot(tt,tmp, '*', 'color', col(j,:));    
end

%% NDF3.9

path2data = '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/data004/data004';
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
    plot(tt+2,tmp, '*', 'color', col(j,:));     
end

%% NDF3.0

path2data = '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/data005/data005';
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
    plot(tt+4,tmp, '*', 'color', col(j,:));     
end



%% NDF2.3

path2data = '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/data007/data007';
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
    plot(tt+6,tmp, '*', 'color', col(j,:));     
end



%% NDF2

path2data = '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/data009/data009';
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
    plot(tt+8,tmp, '*', 'color', col(j,:));     
end


%% NDF1.6

path2data = '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/data011/data011';
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
    plot(tt+10,tmp, '*', 'color', col(j,:));     
end

%% NDF1.3

path2data = '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/data013/data013';
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
    plot(tt+12,tmp, '*', 'color', col(j,:));     
end

%% NDF0.6

path2data = '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/data015/data015';
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
    plot(tt+14,tmp, '*', 'color', col(j,:));     
end

%%
legend('ON parasol', 'OFF parasol', 'ON midget', 'OFF midget')

set(gca, 'xtick', 0.5:2:15, 'xticklabel', {'NDF5.1', 'NDF3.9', 'NDF3', 'NDF2.3', 'NDF2', 'NDF1.6', 'NDF1.3', 'NDF0.6'})

title('2007-02-06-4')