sta1run = load_data('/Volumes/Analysis/2015-10-06-2/d00-13-norefit/data001/data001');
sta1run = load_params(sta1run,'verbose',1);
sta1run = load_sta(sta1run);

sta2run1 = load_data('/Volumes/Analysis/2015-10-06-2/d00-13-norefit/data004/data004');
sta2run1 = load_params(sta2run1,'verbose',1);
sta2run1 = load_sta(sta2run1);

sta2run2 = load_data('/Volumes/Analysis/2015-10-06-2/d00-13-norefit/data011/data011');
sta2run2 = load_params(sta2run2,'verbose',1);
sta2run2 = load_sta(sta2run2);


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

% plot
for cellID = 1:8%128 % 1202
    cellID
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
    
    saveas(gcf, ['/Users/alexth/Desktop/2015-10-06-2_JITTER/svg/Cell_', int2str(cellID), '_datarunID_', int2str(datarunID), '.svg'])
    saveas(gcf, ['/Users/alexth/Desktop/2015-10-06-2_JITTER/tiff/svgCell_', int2str(cellID), '_datarunID_', int2str(datarunID), '.tiff'])
    
    close all
  
end






%% additional stuff
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
figure
colormap gray
imagesc(tmp(:,:,1))

figure
colormap gray
imagesc(inputs_rect(:,:,1))

figure
colormap gray
imagesc(reshape(inputs(:,1), 160, 160))
