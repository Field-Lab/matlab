
datarun=load_data('/Volumes/Analysis/2015-10-29-2/data048-data049/data049/data049');
datarun = load_neurons(datarun);
cell = 3364;
prepped_data = interleaved_data_prep(datarun, 2400, 30, 'cell_spec', cell, 'stimulus_name', '/Volumes/Lab/Users/Nora/new_stim_nora/NSbrownian_code/newrawmovie/gain_control.rawMovie');
cid = get_cell_indices(datarun, cell);


%%
all_spikes = cell2mat(prepped_data.testspikes(:));
spike_frame = floor(all_spikes*120);
spikes = zeros(2400, 1);
for i = 1:2400
    spikes(i) = sum(spike_frame == i); 
end
%spikes = conv(spikes, gausswin(5), 'same');
figure; plot(spikes)
%%
idx = 1:120;
figure; hold on
peak_firing = zeros(10);
for i = 1:10
    image_repeat = spikes(idx);
    plot(image_repeat);
    idx = idx+240;
    peak_firing(i) = max(image_repeat); 
end

%%
% Calculate "expected peak firing"
datarun_class = load_data('/Volumes/Analysis/2015-10-29-2/data048-data049/data048/data048');
datarun_class = load_neurons(datarun_class);
datarun_class = load_sta(datarun_class);

%%
STA = datarun_class.stas.stas{cid};
a = significant_stixels(STA);
STA(~repmat(full(a),[1,1, 3, 30])) = 0.00001;
datarun_class = load_params(datarun_class);
[~, BW_STA] = RGB_to_BW(datarun_class, cid, 'color_movie', STA);
 STA_lin = reshape(BW_STA, [40*80,30]);
[u,s,v] = svd(STA_lin);
time = v(:,1);
space = reshape(u(:,1), [40 80]);
imagesc(space)

%%
% space = imresize(space,4,'nearest');
% BW_STA = permute(BW_STA,[3,2,1]);
% BW_STA = flip(BW_STA, 1);
% BW_STA = flip(BW_STA, 2);
% BW_STA = flip(BW_STA, 3);


%%
% tic
% testmovie = zeros([2400, 80, 40]);
% idx = 1:120;
% for i = 1:10
% testmovie(idx,:,:) =1;
% idx = idx+240;
% end
% testmovie = permute(testmovie, [3, 2, 1]);
% 
% %%
% test = convn(testmovie,space, 'valid');
% drive = conv(squeeze(test), time, 'full');
% drive = drive(1:2400);

%%
testmovie = permute(prepped_data.testmovie, [3,2,1]);
testmovie_down = zeros(40, 80, 2400);
idx_i = 1:2;
for i = 1:40
    idx_i=idx_i+2;
    idx_j = 1:2;
    for j=1:80
       idx_j = idx_j+2;
       testmovie_down(i,j,:) = sum(sum(prepped_data.testmovie(:,idx_j, idx_i),3),2);
    end
end

%%
linear_drive = convn(testmovie_down, -space, 'valid');
linear_drive = conv(squeeze(linear_drive), time, 'full');
linear_drive = linear_drive(1:2400);
plot(linear_drive)

%%
%linear_drive = [zeros(29,1); linear_drive];
plot(10*(exp(linear_drive/10^5)-1));
hold on; plot(spikes/max(spikes));
for i = 1:20; hold on; plot(i*[120,120], [0 1], 'yellow'); end

%%
figure;
hold on
idx = 1:120;
for i = 1:10
    image_repeat = linear_drive(idx);
    plot(image_repeat+i);
    idx = idx+240;
    %peak_firing(i) = max(image_repeat); 
end
