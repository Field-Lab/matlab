
datarun = load_data('/Volumes/Analysis/2015-10-06-2/d00-13-norefit/data008/data008');
datarun = load_params(datarun,'verbose',1);
% datarun = load_sta(datarun);
datarun = load_neurons(datarun);

triggers = datarun.triggers;
movie_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-2-8-0.48-11111-160x160.xml';



datarun = load_data('/Volumes/Acquisition/Analysis/2016-02-17-7/data003/data003');
datarun = load_params(datarun,'verbose',1);
% datarun = load_sta(datarun);
datarun = load_neurons(datarun);
triggers = datarun.triggers;

movie_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-2-6-0.48-11111-400x300-60.35.xml';
tic
[mov,height,width,duration,refresh] = get_movie(movie_file,triggers, Inf);
toc

rgb_flag  = 'bw'; bin_flag = 'bin'; rgb = 0.48; seed = 11111; width = 400; height = 300; probability = 1; frames = 20734;
tic
my_movie = compute_raw_wn_movie(rgb_flag, bin_flag, rgb, seed, width, height, probability, frames);
toc

a= my_movie(:,:,1,1);
b = mov(:,:,1,1);

figure
imagesc(a');
axis([0 20 0 20])
figure
imagesc(b)
axis([0 20 0 20])

rgb_flag  = 'bw'; bin_flag = 'bin'; rgb = 0.48; seed = 11111; width = 2; height = 2; probability = 1; frames = 27026;
tic
my_movie = compute_raw_wn_movie(rgb_flag, bin_flag, rgb, seed, width, height, probability, frames);
toc




[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-2-8-0.48-11111-160x160.xml');

parameters.seed = 11111;
stimulus.stixel_width = 2;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = 160;  
parameters.field_height = 160;

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
clear full_inputs inputs

tmp = (full_inputs1+0.48)/0.96;
tmp = uint8(tmp);
save('/Users/alexth/Desktop/2015-10-06-2_JITTER/inputs_jitter.mat', 'tmp', '-v7.3');

%% main stuff
sta_params.length = 5;
sta_params.offset = 0;
% all_sta = cell(1,128);
load('/Users/alexth/Desktop/2015-10-06-2_JITTER/correct_jitter_sta.mat');

% calc
for cellID = 34:128 
    cellID
    spikes = datarun.spikes{cellID};
    spikes=ceil((spikes-datarun.triggers(1))*1000/(refresh)); % spikes in frames
    
    spikes(spikes>size(full_inputs1,2)) = [];
    spikes(spikes<sta_params.length-sta_params.offset) = [];
    
    
    % jittered inputs
    spikes_tmp = spikes;
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
    sta = reshape(sta, 322,322,sta_params.length);
    all_sta{cellID} = sta;
    
    save('/Users/alexth/Desktop/2015-10-06-2_JITTER/correct_jitter_sta.mat', 'all_sta');

    datarunID = datarun.cell_ids(cellID);
    figure
    set(gcf, 'position', [-1680          97        1032        1008])
    
    sta1 = squeeze(sta1run.stas.stas{cellID});
    sta1 = sta1(:,:,4);
    
    sta2 = squeeze(sta2run1.stas.stas{cellID});
    sta2 = imresize(sta2(:,:,4),2,'method', 'nearest');
    
    sta3 = squeeze(sta2run2.stas.stas{cellID});
    sta3 = imresize(sta3(:,:,4),2,'method', 'nearest');
    
    sta_jitter = all_sta{cellID}(:,:,2);
    
    mean_sd = mean(sta1run.stas.fits{cellID}.sd);
    centers = sta1run.stas.fits{cellID}.mean;
    
    subplot(2,2,1)
    colormap gray
    imagesc(sta1)
    axis([centers(1)-mean_sd*3 centers(1)+mean_sd*3 centers(2)-mean_sd*3 centers(2)+mean_sd*3])
    title(['cell ', int2str(cellID), ', datarun ID ', int2str(datarunID), ', BW-1-8, data001'], 'Interpreter', 'none')
    set(gca, 'DataAspectRatio', [1 1 1])
    
    subplot(2,2,2)
    colormap gray
    imagesc(sta2)
    axis([centers(1)-mean_sd*3 centers(1)+mean_sd*3 centers(2)-mean_sd*3 centers(2)+mean_sd*3])
    title(['cell ', int2str(cellID), ', datarun ID ', int2str(datarunID), ', BW-2-8, data004'], 'Interpreter', 'none')
    set(gca, 'DataAspectRatio', [1 1 1])
    
    subplot(2,2,4)
    colormap gray
    imagesc(sta3)
    axis([centers(1)-mean_sd*3 centers(1)+mean_sd*3 centers(2)-mean_sd*3 centers(2)+mean_sd*3])
    title(['cell ', int2str(cellID), ', datarun ID ', int2str(datarunID), ', BW-2-8, data011'], 'Interpreter', 'none')
    set(gca, 'DataAspectRatio', [1 1 1])
    
    subplot(2,2,3)
    colormap gray
    imagesc(sta_jitter)
    axis([centers(1)-mean_sd*3 centers(1)+mean_sd*3 centers(2)-mean_sd*3 centers(2)+mean_sd*3])
    title(['cell ', int2str(cellID), ', datarun ID ', int2str(datarunID), ', BW-2-8 JITTER, data008'], 'Interpreter', 'none')
    set(gca, 'DataAspectRatio', [1 1 1])
    
    drawnow
    
    saveas(gcf, ['/Users/alexth/Desktop/2015-10-06-2_JITTER/svg/Cell_', int2str(cellID), '_datarunID_', int2str(datarunID), '.svg'])
    saveas(gcf, ['/Users/alexth/Desktop/2015-10-06-2_JITTER/tiff/svgCell_', int2str(cellID), '_datarunID_', int2str(datarunID), '.tiff'])
    
    close all
 

end


save('/Users/alexth/Desktop/correct_jitter_sta.mat', 'all_sta');
save('/Users/alexth/Desktop/inputs_jitter.mat', 'full_inputs1', '-v7.3');
