datarun = load_data('/Volumes/Analysis/2015-10-06-2/d00-13-norefit/data008/data008');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = load_neurons(datarun);

%%BW-2-8-0.48-11111 JITTER 320x320 field LISP
[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-2-8-0.48-11111-160x160.xml');

parameters.seed = 11111;
stimulus.stixel_width = 2;
stimulus.stixel_height = 2;

stimulus.rng_init.state = parameters.seed;
stimulus.rng_init.seed = parameters.seed;
stimulus.rng_init.state = Init_RNG_JavaStyle(stimulus.rng_init.seed);
stimulus.jitter.state = stimulus.rng_init.state;

full_inputs = zeros(322,322,size(inputs,2));
for i=1:size(inputs,2)
    if mod(i,1000)==0
        i
    end
    tmp = reshape(inputs(:,i), 160, 160);    
    tmp = imresize(tmp,2, 'method', 'nearest');
    jitterX = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_width) - stimulus.stixel_width/2;
    jitterY = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_height) - stimulus.stixel_height/2;
    
    full_inputs(2+jitterX:321+jitterX,2+jitterY:321+jitterY,i) = tmp;
end

full_inputs1 = reshape(full_inputs,322*322,size(full_inputs,3));



cellID = 14; % 1202


sta_params.length = 5;
sta_params.offset = 0;

spikes = datarun.spikes{cellID};
spikes=ceil((spikes-datarun.triggers(1))*1000/(refresh)); % spikes in frames

spikes(spikes>size(inputs,2)) = []; 
spikes(spikes<sta_params.length-sta_params.offset) = [];
% get spike rate
spikes_tmp = spikes;
spike_rate=zeros(size(inputs,2),1);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp

% normal inputs
spikes_tmp = spikes; 
sta=zeros(size(inputs,1),sta_params.length);
nspikes = numel(spikes_tmp);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    for j=1:sta_params.length
        sta(:,sta_params.length-j+1) = sta(:,sta_params.length-j+1)...
            + sum(inputs(:,spikes_tmp(ia) - sta_params.length + j + sta_params.offset),2);
    end
    spikes_tmp(ia)=[];
end
sta = sta/nspikes;

figure
for i=1:5
    subplot(2,3,i)
    colormap gray
    imagesc(reshape(sta(:,i), 160, 160))
end

% jittered inputs
spikes_tmp = spikes(end/3:end*2/3); 
sta=zeros(size(full_inputs1,1),sta_params.length);
nspikes = numel(spikes_tmp);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    length(ia)
    for j=1:sta_params.length
        sta(:,sta_params.length-j+1) = sta(:,sta_params.length-j+1)...
            + sum(full_inputs1(:,spikes_tmp(ia) - sta_params.length + j + sta_params.offset),2);
    end
    spikes_tmp(ia)=[];
end
sta = sta/nspikes;

figure
for i=1:5
    subplot(2,3,i)
    colormap gray
    imagesc(reshape(sta(:,i), 322, 322))
end



figure
colormap gray
imagesc(tmp(:,:,1))

figure
colormap gray
imagesc(inputs_rect(:,:,1))

figure
colormap gray
imagesc(reshape(inputs(:,1), 160, 160))
