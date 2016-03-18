local_path = '/Volumes/Analysis/';

%% 2016-02-17-4
%SC run 1
datarun = load_data([local_path, '2016-02-17-4/d00-05-norefit/data001/data001']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = load_neurons(datarun);

% collect STAs
all_sta = zeros(300,400,length(datarun.cell_ids));
cnt = 1;
pols = zeros(1,length(datarun.cell_ids));
max_pics = pols;
for i = 1:length(datarun.cell_ids)
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
    all_sta(:,:,cnt)=double(squeeze(sta));
    tmp = sta*pols(cnt);
    max_pics(cnt) = max(tmp(:));
    cnt = cnt+1;
end

% noise maps
cell_inds = get_cell_indices(datarun,{3,4,5});
noise_sta = mean(all_sta(:,:,cell_inds),3);
load('/Volumes/Analysis/2016-02-17-4/noise_BW_2_6_11111.mat')
noise_stimul = p;

denoised_sta = zeros(300,400,6,length(datarun.cell_ids));
for i = 1:length(datarun.cell_ids)
    denoised_sta(:,:,:,i) = squeeze(datarun.stas.stas{i})-repmat(noise_stimul,1,1,6);
end
denoised_sta = noise_sta;
save('/Volumes/Analysis/2016-03-17-0/data002_denoised_sta', 'denoised_sta')

% mean STA of several cells
cell_inds = get_cell_indices(datarun,{2});
tmp = mean(all_sta(:,:,cell_inds(1:5)),3);
figure
set(gcf, 'position', [-1753 529 1679 536])
ax1 = subplot(1,3,1);
imagesc(tmp);
set(ax1, 'DataAspectRatio', [1 1 1])
title(['raw mean, std ',num2str(robust_std(tmp(:))) ])
ax2 = subplot(1,3,2);
tmp2 = tmp- noise_sta;
imagesc(tmp2);
set(ax2, 'DataAspectRatio', [1 1 1])
title(['denoised mean cells, std ',num2str(robust_std(tmp2(:)))])
ax3 = subplot(1,3,3);
tmp2 = tmp- noise_stimul;
imagesc(tmp2);
set(ax3, 'DataAspectRatio', [1 1 1])
title(['denoised mean stimulus, std ',num2str(robust_std(tmp2(:)))])
linkaxes([ax1, ax2, ax3], 'xy')



cell_inds = get_cell_indices(datarun,{2});
for i=1:10
    tmp = sum(all_sta(:,:,cell_inds(i)),3);
    tmp1 = tmp-noise_sta;
    tmp2 = tmp - noise_stimul;
    figure
    set(gcf, 'position', [-1870         550        1871         521])
    ax1 = subplot(1,3,1);
    colormap gray
    imagesc(tmp1)
    title(['denoised with cells , ', num2str(robust_mean(tmp1(:))), ' + ', num2str(robust_std(tmp1(:)))])
    ax2 = subplot(1,3,2);
    colormap gray
    imagesc(tmp2)
    title(['denoised with stimulus, ', num2str(robust_mean(tmp2(:))), ' + ', num2str(robust_std(tmp2(:)))])
    ax3 = subplot(1,3,3);
    colormap gray
    imagesc(tmp)
    title(['raw, ', num2str(robust_mean(tmp(:))), ' + ', num2str(robust_std(tmp(:)))])
    linkaxes([ax1,ax2,ax3],'xy')
end


%% sig stixels approach

threshold = 8;
% noise maps
cell_inds = get_cell_indices(datarun,{1,2,3,4});
noise_sta = 0;
for i=1:length(cell_inds)
    tmp = all_sta(:,:,cell_inds(i));
    p = robust_std(tmp(:));
    [a,b] = find(abs(tmp)>p*threshold);
    for j=1:length(a)
        tmp(a(j), b(j)) = 0;
    end
    noise_sta = noise_sta+tmp;
end
noise_sta = noise_sta/length(cell_inds);

% mean STA of several cells
cell_inds = get_cell_indices(datarun,{4});
tmp = mean(all_sta(:,:,cell_inds(1:end)),3);
figure
set(gcf, 'position', [-1753 529 1679 536])
ax1 = subplot(1,2,1);
imagesc(tmp);
set(ax1, 'DataAspectRatio', [1 1 1])
title(['raw mean, std ',num2str(robust_std(tmp(:))) ])
ax2 = subplot(1,2,2);
tmp2 = tmp-noise_sta;
imagesc(tmp2);
set(ax2, 'DataAspectRatio', [1 1 1])
title(['denoised mean, std ',num2str(robust_std(tmp2(:)))])
linkaxes([ax1, ax2], 'xy')


%% 2016-02-17-7
%SC run 1
datarun1 = load_data([local_path, '2016-02-17-7/d01-13-norefit/data002/data002']);
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = load_neurons(datarun1);

% collect STAs
all_sta1 = zeros(300,400,length(datarun1.cell_ids));
cnt = 1;
pols = zeros(1,length(datarun1.cell_ids));
max_pics = pols;
for i = 1:length(datarun1.cell_ids)
    tmp = datarun1.stas.stas{i};
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
    all_sta1(:,:,cnt)=double(squeeze(sta));
    tmp = sta*pols(cnt);
    max_pics(cnt) = max(tmp(:));
    cnt = cnt+1;
end



% noise maps
cell_inds = get_cell_indices(datarun1,{2,3,4});
noise_sta = mean(all_sta1(:,:,cell_inds),3);

load('/Volumes/Analysis/2016-02-17-4/noise_BW_2_6_11111.mat')
noise_stimul = r;

threshold = 8;
noise_sta_sigstix = 0;
for i=1:length(cell_inds)
    tmp = all_sta1(:,:,cell_inds(i));
    p = robust_std(tmp(:));
    [a,b] = find(abs(tmp)>p*threshold);
    for j=1:length(a)
        tmp(a(j), b(j)) = 0;
    end
    noise_sta_sigstix = noise_sta_sigstix+tmp;
end
noise_sta_sigstix = noise_sta_sigstix/length(cell_inds);


% mean STA of several cells
cell_inds = get_cell_indices(datarun1,{4});
tmp = mean(all_sta1(:,:,cell_inds(1:5)),3);
figure
set(gcf, 'position', [-1753 529 1679 536])
ax1 = subplot(1,3,1);
imagesc(tmp);
set(ax1, 'DataAspectRatio', [1 1 1])
title(['raw mean, std ',num2str(robust_std(tmp(:))) ])
ax2 = subplot(1,3,2);
tmp2 = tmp- noise_sta;
imagesc(tmp2);
set(ax2, 'DataAspectRatio', [1 1 1])
title(['denoised mean with cells, std ',num2str(robust_std(tmp2(:)))])
ax3 = subplot(1,3,3);
tmp2 = tmp- noise_stimul;
imagesc(tmp2);
set(ax3, 'DataAspectRatio', [1 1 1])
title(['denoised mean with stimul, std ',num2str(robust_std(tmp2(:)))])
linkaxes([ax1, ax2, ax3], 'xy')



cell_inds = get_cell_indices(datarun1,{2});
for i=1:10
    tmp = sum(all_sta1(:,:,cell_inds(i)),3);
    tmp1 = tmp-noise_sta;
    tmp2 = tmp - noise_stimul;
    figure
    set(gcf, 'position', [-1870         550        1871         521])
    ax1 = subplot(1,3,1);
    colormap gray
    imagesc(tmp1)
    title(['denoised with cells , ', num2str(robust_mean(tmp1(:))), ' + ', num2str(robust_std(tmp1(:)))])
    ax2 = subplot(1,3,2);
    colormap gray
    imagesc(tmp2)
    title(['denoised with stimulus, ', num2str(robust_mean(tmp2(:))), ' + ', num2str(robust_std(tmp2(:)))])
    ax3 = subplot(1,3,3);
    colormap gray
    imagesc(tmp)
    title(['raw, ', num2str(robust_mean(tmp(:))), ' + ', num2str(robust_std(tmp(:)))])
    linkaxes([ax1,ax2,ax3],'xy')
end



%%

p = mean(all_sta([1:10 291:300],[1:10 391:400],get_cell_indices(datarun, {1,2,3,4})), 3);
k = mean(all_sta1([1:10 291:300],[1:10 391:400],get_cell_indices(datarun1, {2,3,4})), 3);
figure
plot(p(:))
hold on
plot(k(:))
title(['Correlation p = ', num2str(corr(p(:), k(:)))])

cell_inds = get_cell_indices(datarun, {1,2,3,4});
all_spikes = round(sort(cell2mat(datarun.spikes(cell_inds)))*1000/99.4200);
all_spikes(all_spikes==0) = [];
spike_rate = zeros(max(all_spikes),1);
while ~isempty(all_spikes)
    [~, ia, ~] = unique(all_spikes);
    spike_rate(all_spikes(ia))=spike_rate(all_spikes(ia))+1;
    all_spikes(ia)=[];
end

cell_inds1 = get_cell_indices(datarun1, {2,3,4});
all_spikes1 = round(sort(cell2mat(datarun1.spikes(cell_inds1)))*1000/99.4200);
all_spikes1(all_spikes1==0) = [];
spike_rate1 = zeros(max(all_spikes1),1);
while ~isempty(all_spikes1)
    [~, ia, ~] = unique(all_spikes1);
    spike_rate1(all_spikes1(ia))=spike_rate1(all_spikes1(ia))+1;
    all_spikes1(ia)=[];
end

figure
plot(spike_rate)
hold on
plot(spike_rate1)
title(['Correlation p = ', num2str(corr(spike_rate(1:length(spike_rate1)), spike_rate1))])




p = mean(all_sta([1:10 291:300],[1:10 391:400],get_cell_indices(datarun, {3})), 3);
k = mean(all_sta([1:10 291:300],[1:10 391:400],get_cell_indices(datarun1, {4})), 3);
figure
plot(p(:))
hold on
plot(k(:))


p = rand(100,1);
q = rand(100,1);

corr(p,q)

