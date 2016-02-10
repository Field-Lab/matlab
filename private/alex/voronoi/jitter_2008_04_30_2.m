

datarun = load_data('/Volumes/Analysis/2008-04-30-2/data007/data007');
datarun = load_params(datarun,'verbose',1);
% datarun = load_sta(datarun);
datarun = load_neurons(datarun);

sta_params.length = 5;
sta_params.offset = 0;
all_sta = cell(1,length(datarun.cell_ids));

%%BW-2-8-0.48-11111 JITTER 320x320 field LISP, stixel field 128x64
[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-5-6-0.48-11111.xml');

spike_array = uint8(zeros(length(datarun.cell_ids), duration));
for cellID = 1:length(datarun.cell_ids)    
    spikes = datarun.spikes{cellID};
    spikes=ceil((spikes-datarun.triggers(1))*1000/(refresh)); % spikes in frames
    spikes(spikes>duration) = [];
    spikes(spikes<sta_params.length-sta_params.offset) = [];    
    spike_array(cellID, spikes) = 1;
end

load('/Volumes/Analysis/2008-04-30-2/jitter/shifts')
% parameters.seed = 11111;
% stimulus.stixel_width = 5;
% stimulus.stixel_height = 5;

% stimulus.rng_init.state = parameters.seed;
% stimulus.rng_init.seed = parameters.seed;
% stimulus.rng_init.state = Init_RNG_JavaStyle(stimulus.rng_init.seed);
% stimulus.jitter.state = stimulus.rng_init.state;

full_inputs = zeros(640+4,320+4,sta_params.length);
for i=1:sta_params.length-1
    tmp = reshape(inputs(:,i), 128, 64);
    tmp = imresize(tmp,5, 'method', 'nearest');
    jitterX = shifts(1,i);
    jitterY = shifts(2,i);
%     jitterX = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_width) - floor(stimulus.stixel_width/2);
%     jitterY = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_height) - floor(stimulus.stixel_height/2);
    full_inputs(3+jitterX:640+2+jitterX,3+jitterY:322+jitterY,1+i) = tmp;    
end

sta = zeros(644,324, sta_params.length,length(datarun.cell_ids) );
for i=sta_params.length:duration 
    i
    tmp = reshape(inputs(:,i), 128, 64);
    tmp = imresize(tmp,5, 'method', 'nearest');
    jitterX = shifts(1,i);
    jitterY = shifts(2,i);
%     jitterX = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_width) - floor(stimulus.stixel_width/2);
%     jitterY = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_height) - floor(stimulus.stixel_height/2);
    
    full_inputs = circshift(full_inputs,-1,3);
    full_inputs(3+jitterX:640+2+jitterX,3+jitterY:322+jitterY,sta_params.length) = tmp;   
    
    a = find(spike_array(:,i));
    sta(:,:,:,a) =  sta(:,:,:,a) + repmat(full_inputs, 1, 1, 1, length(a));    
end

for i = 1:length(datarun.cell_ids)
    sta(:,:,:,i) =  sta(:,:,:,i) / nnz(spike_array(i,:));
end

save('/Volumes/Analysis/2008-04-30-2/jitter/correct_jitter_sta.mat', 'sta');
