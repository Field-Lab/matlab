
col = [1 0 0; 0 0 1; 0.7 0 0.7; 0.1 0.5 0.1];
% 519 system (320x320)

% data005 (5, BW-20-10), data008 (4, BW-10-8)
% data011 (3, BW-16-8), data014 (3, BW-8-8)
% data015 (2, BW-16-6), data018 (2, BW-8-6)
% data019 (1, RGB-16-6), data022 (1, BW-8-6)
% data023 (0, RGB-16-4), data026 (0, BW-4-2)
% data027 (0, RGB-4-2)

chrono_order = [0:2:22];
% nd_order = [14 12 10 4 0 16 6 2 8 18];

ord = chrono_order;

%% NDF5
path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data005/data005';
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

%% NDF4

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data008/data008';
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

%% NDF3

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data011/data011';
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



%% NDF3

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data014/data014';
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



%% NDF2

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data015/data015';
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


%% NDF2

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data018/data018';
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

%% NDF1

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data019/data019';
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

%% NDF1

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data022/data022';
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


%% NDF0

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data023/data023';
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

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data026/data026';
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

%% NDF0

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data027/data027';
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
    plot(tt+ord(11),tmp, '*', 'color', col(j,:));     
end



%%
legend('ON parasol', 'OFF parasol', 'ON midget', 'OFF midget')

% chrono
set(gca, 'xtick', 0.5:2:22, 'xticklabel', {'NDF5\newlineBW-20', 'NDF4\newlineBW-10',...
    'NDF3\newlineBW-16', 'NDF3\newlineBW-8', 'NDF2\newlineBW-16', 'ND2\newlineBW-8',...
    'NDF1\newlineRGB-16', 'NDF1\newlineBW-8', 'NDF0\newlineRGB-16', 'NDF0\newlineBW-4',...
    'NDF0\newlineRGB-4'})
title('2015-03-09-2, chrono order')

% % ndf
% set(gca, 'xtick', 0.5:2:19, 'xticklabel', {'NDF3.3', 'NDF2.6', 'NDF1.6', 'NDF1.3', 'NDF1.3', 'NDF1.0', 'NDF0.6', 'NDF0', 'NDF0', 'NDF0'})
% title('2007-02-06-4, NDF order')


%%
%%%%%%%%%%%%%%%%%%%%%%%% Only 16 stixel %%%%%%%%%%%%%%%%%%%%%%%%

%%

col = [1 0 0; 0 0 1; 0.7 0 0.7; 0.1 0.5 0.1];
% 519 system (320x320)

% data011 (3, BW-16-8),
% data015 (2, BW-16-6),
% data019 (1, RGB-16-6),
% data023 (0, RGB-16-4),

chrono_order = [0:2:22];
% nd_order = [14 12 10 4 0 16 6 2 8 18];
ord = chrono_order;



%% NDF3

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data011/data011';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

figure
a = zeros(100,4);
for j=1:4
    cellInds = get_cell_indices(datarun, {j});
    for cellID=cellInds
        a(cellID,j) = fit_area(datarun, cellID);
    end
    tmp = a(a(:,j)>0,j);
    tt = rand(size(tmp));
    plot(tt+ord(1),tmp, '*', 'color', col(j,:));
    hold on
end



%% NDF2

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data015/data015';
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


%% NDF1

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data019/data019';
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


%% NDF0

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data023/data023';
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




%%
legend('ON parasol', 'OFF parasol', 'ON midget', 'OFF midget')

% chrono
set(gca, 'xtick', 0.5:2:22, 'xticklabel', {'NDF3\newlineBW-16', 'NDF2\newlineBW-16',...
    'NDF1\newlineRGB-16', 'NDF0\newlineRGB-16'})
title('2015-03-09-2, only stixel 16')
