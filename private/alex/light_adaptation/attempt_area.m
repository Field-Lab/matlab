%% NDF3

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data011/data011';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

c = datarun.stas.stas{get_cell_indices(datarun, 263)};

c = double(squeeze(c));

figure
colormap gray
imagesc(c(:,:,27))

figure
plot(datarun.vision.timecourses(get_cell_indices(datarun, 1921)).g )


%% NDF3

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data014/data014';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

a = datarun.stas.stas{get_cell_indices(datarun, 263)};

a = double(squeeze(a));

figure
colormap gray
imagesc(a(:,:,27))

figure
plot(datarun.vision.timecourses(get_cell_indices(datarun, 1921)).g )

%% NDF2
clear datarun
path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data018/data018';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

b = datarun.stas.stas{get_cell_indices(datarun, 263)};

b= double(squeeze(b));

figure
colormap gray
imagesc(b(:,:,27))

figure
plot(datarun.vision.timecourses(get_cell_indices(datarun, 1921)).g )


%% NDF1
clear datarun
path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data022/data022';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

c = datarun.stas.stas{get_cell_indices(datarun, 263)};

c= double(squeeze(c));

figure
colormap gray
imagesc(c(:,:,28))

figure
plot(datarun.vision.timecourses(get_cell_indices(datarun, 1921)).g )

%%
tmp = a(:,:,27)/max(a(:));
tmp1 = b(:,:,27)/max(b(:));
tmp2 = c(:,:,28)/max(c(:));

new_map = tmp-tmp1;


figure
colormap gray
imagesc(new_map)

k = zeros(40,40,3);
% k(:,:,1) = tmp;
k(:,:,2) = tmp;
k(:,:,3) = tmp1;
figure
imagesc(k)


%% NDF2

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data015/data015';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);


a = datarun.stas.stas{get_cell_indices(datarun, 263)};

a = double(squeeze(a));


clear datarun
path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data018/data018';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all');
datarun = get_sta_fits_from_vision(datarun);

b = datarun.stas.stas{get_cell_indices(datarun, 263)};

b= double(squeeze(b));

tmp = a(:,:,27)/max(a(:));
tmp1 = b(:,:,27)/max(b(:));

tmp = imresize(tmp,2);

k = zeros(40,40,3);
% k(:,:,1) = tmp;
k(:,:,2) = tmp;
k(:,:,3) = tmp1;
figure
imagesc(k)



c1 = imresize(c(:,:,27), 2);
figure
colormap gray
imagesc(c1)

a1 = a(:,:,27);
b1 = b(:,:,27);

t = zeros(40,40,3);
t(:,:,1) = a1/max(a1(:));
t(:,:,2) = c1/max(c1(:));
figure
subplot(2,2,1)
imagesc(t)
title('NDF3 size 8(red), NDF3 size 16(green)')

t = zeros(40,40,3);
t(:,:,1) = a1/max(a1(:));
t(:,:,2) = b1/max(b1(:));
subplot(2,2,2)
imagesc(t)
title('NDF3 size 8(red), NDF2 size 8(green)')

t = zeros(40,40,3);
t(:,:,1) = c1/max(c1(:));
t(:,:,2) = b1/max(b1(:));
subplot(2,2,3)
imagesc(t)
title('NDF3 size 16(red), NDF2 size 8(green)')



% nearest



c1 = imresize(c(:,:,27), 2, 'method', 'nearest');
figure
colormap gray
imagesc(c1)

a1 = a(:,:,27);
b1 = b(:,:,27);

t = zeros(40,40,3);
t(:,:,1) = a1/max(a1(:));
t(:,:,2) = c1/max(c1(:));
figure
subplot(2,2,1)
imagesc(t)
title('NDF3 size 8(red), NDF3 size 16(green)')

t = zeros(40,40,3);
t(:,:,1) = a1/max(a1(:));
t(:,:,2) = b1/max(b1(:));
subplot(2,2,2)
imagesc(t)
title('NDF3 size 8(red), NDF2 size 8(green)')

t = zeros(40,40,3);
t(:,:,1) = c1/max(c1(:));
t(:,:,2) = b1/max(b1(:));
subplot(2,2,3)
imagesc(t)
title('NDF3 size 16(red), NDF2 size 8(green)')

