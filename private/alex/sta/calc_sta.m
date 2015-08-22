cd('/Users/alexth/test4/matlab/private/alex/sta/')
mex sta.c

init_state = Init_RNG_JavaStyle(11111);
width = 10;
height = 10;
map = [];
map_back_rgb = uint8([128 128 128]);

noise_type = 0;
n_bits = 1;
probability = 1;
sta_length = 5;

spikes = {[10, 30], [15]};
frames = max(cellfun(@max,spikes))+sta_length+1;
% spikes = {[2 4 6], [8]};
ncells = size(spikes,2);
bin_spikes_array = zeros(ncells, frames+1);
for i=1:ncells
    bin_spikes_array(i,spikes{i}) = bin_spikes_array(i,spikes{i})+1; % NB: if 2 spikes per frame?
end
bin_spikes_array = uint8(bin_spikes_array)';

rgb_vect = [1 1 1]*0.48;
back_rgb = [1 1 1]*0.5;
tmp = [1 1 1; -1 -1 -1] .* repmat(rgb_vect,2,1)  + repmat(back_rgb,2,1);
tmp = uint8(round(255 * tmp))';
lut = tmp(:);




my_sta = sta(init_state, width, height, lut, map, map_back_rgb,noise_type,n_bits,probability, frames, ncells, sta_length, bin_spikes_array);
size(my_sta)

my_sta = squeeze(my_sta(:,:,1)); 
size(my_sta)
my_sta = reshape(my_sta,3,width,height, sta_length); % colors, width, height, frames
nnz(squeeze(my_sta(1,:,:,:))-squeeze(my_sta(3,:,:,:)))



%% real cells

%data011 BW-16-8-0.48-11111 

datarun = load_data('/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data011-from-d05-d27/data011-from-d05-d27');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarunID=find(datarun.cell_ids==47);

% convert to frames
[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-16-8-0.48-11111-20x20.xml');


width = 20;
height = 20;
map = [];
map_back_rgb = uint8([128 128 128]);
rgb_vect = [1 1 1]*0.48;
back_rgb = [1 1 1]*0.5;
tmp = [1 1 1; -1 -1 -1] .* repmat(rgb_vect,2,1)  + repmat(back_rgb,2,1);
tmp = uint8(round(255 * tmp))';
lut = tmp(:);

noise_type = 0;
n_bits = 1;
probability = 1;
sta_length = 30;

frames = 28000;%max(cellfun(@max,spikes))+sta_length+1;
ncells = length(datarun.cell_ids);
all_stas = zeros(1200,30,ncells);
nspikes = zeros(ncells,1);
tic
parfor i=1:ncells
    bin_spikes_array = zeros(1, frames+1);
    spikes=ceil((datarun.spikes{i}-datarun.triggers(1))*1000/(refresh)); % spikes in frames
    while ~isempty(spikes)
        [~, ia, ~] = unique(spikes);
        bin_spikes_array(spikes(ia))=bin_spikes_array(spikes(ia))+1;
        spikes(ia)=[];
    end
    bin_spikes_array = uint8(bin_spikes_array)';
    nspikes(i) = sum(bin_spikes_array);
    init_state = Init_RNG_JavaStyle(11111);
    all_stas(:,:,i) = sta(init_state, width, height, lut, map, map_back_rgb,noise_type,n_bits,probability, frames, 1, sta_length, bin_spikes_array);

end
toc







frames = 28000;%max(cellfun(@max,spikes))+sta_length+1;
ncells = length(datarun.cell_ids);
all_stas = zeros(1200,30,ncells);
nspikes = zeros(ncells,1);

bin_spikes_array = zeros(ncells, frames+1);
for i=1:ncells
    spikes=ceil((datarun.spikes{i}-datarun.triggers(1))*1000/(refresh)); % spikes in frames        
    while ~isempty(spikes)
        [~, ia, ~] = unique(spikes);
        bin_spikes_array(i,spikes(ia))=bin_spikes_array(i,spikes(ia))+1;
        spikes(ia)=[];
    end    
   
end
bin_spikes_array = uint8(bin_spikes_array)';
nspikes = sum(bin_spikes_array);

tic
parfor i=1:ncells
    init_state = Init_RNG_JavaStyle(11111);
    all_stas(:,:,i) = sta(init_state, width, height, lut, map, map_back_rgb,noise_type,n_bits,probability, frames, 1, sta_length, bin_spikes_array(:,i));
end
toc




cellID = 1;
tmp_sta = squeeze(all_stas(:,:,cellID))/nspikes(cellID); 
size(tmp_sta)
tmp_sta = reshape(tmp_sta,3,width,height, sta_length); % colors, width, height, frames
size(tmp_sta)
nnz(squeeze(tmp_sta(1,:,:,:))-squeeze(tmp_sta(3,:,:,:)))

figure
for i=1:30
    subplot(6,5,i)
    tmp = squeeze(tmp_sta(1,:,:,i));
    colormap gray
    imagesc(tmp)
end



% test for memory
ncells = 50;
bin_spikes_array = zeros(1, frames+1);
for i=1:ncells    
    spikes=ceil((datarun.spikes{i}-datarun.triggers(1))*1000/(refresh)); % spikes in frames
    while ~isempty(spikes)
        %         bin_spikes_array(i,spikes{i}) = bin_spikes_array(i,spikes{i})+1; % NB: if 2 spikes per frame?
        [~, ia, ~] = unique(spikes);
        bin_spikes_array(i,spikes(ia))=bin_spikes_array(spikes(ia))+1;
        spikes(ia)=[];
    end    
end
bin_spikes_array = uint8(bin_spikes_array)';
init_state = Init_RNG_JavaStyle(11111);
tic
my_sta = sta(init_state, width, height, lut, map, map_back_rgb,noise_type,n_bits,probability, frames, ncells, sta_length, bin_spikes_array);
toc




size(my_sta)
cellID = 3;

tmp_sta = squeeze(my_sta(:,:,cellID))/sum(bin_spikes_array(:,cellID)); 
size(tmp_sta)
tmp_sta = reshape(tmp_sta,3,width,height, sta_length); % colors, width, height, frames
size(tmp_sta)
nnz(squeeze(tmp_sta(1,:,:,:))-squeeze(tmp_sta(3,:,:,:)))

figure
for i=1:30
    subplot(6,5,i)
    tmp = squeeze(tmp_sta(1,:,:,i));
    colormap gray
    imagesc(tmp)
end
figure
tmp=squeeze(tmp_sta(1,15,8,end:-1:1));
plot(tmp)
hold on
tmp=squeeze(tmp_sta(1,14,8,end:-1:1));
plot(tmp)


figure

    tmp = squeeze(my_sta(1,:,:,3));
    colormap gray
    imagesc(tmp)




