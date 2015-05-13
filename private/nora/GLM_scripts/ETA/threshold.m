function [PA, NA] = threshold(xval, movie)

n_blocks = length(movie);
stim_size = size(movie{1,1}.matrix);

n_spikes_pos = 0;
n_spikes_neg = 0;
total_spikes = 0;

PA = zeros(stim_size(1), stim_size(2), 30);
NA = zeros(stim_size(1), stim_size(2), 30);

for i_block = 1:n_blocks
    mov = movie{i_block}.matrix;
    prob_spike = 1 - exp(-xval{i_block}.glm_rateperbin);
    thresh_pos_spikes = find(logical(xval{i_block}.rasters.recorded) & (prob_spike < 0.7) & (prob_spike > 0.3));
    thresh_neg_spikes = find(~logical(xval{i_block}.rasters.recorded) & (prob_spike < 0.7) & (prob_spike > 0.3));
    total_spikes = total_spikes + sum((prob_spike < 0.7) & (prob_spike > 0.3));
    n_spikes_pos = n_spikes_pos + length(thresh_pos_spikes);
    n_spikes_neg = n_spikes_neg + length(thresh_neg_spikes);
    
    for i_spike = 1:length(thresh_pos_spikes)
        frame = ceil(thresh_pos_spikes(i_spike)/10);
        if frame>29
            PA = PA + double(mov(:,:,(frame-29):frame));
        end
    end
    for i_spike = 1:length(thresh_neg_spikes)
        frame = ceil(thresh_neg_spikes(i_spike)/10);
        if frame>29
            NA = NA + double(mov(:,:,(frame-29):frame));
        end
    end
    
end

min1 = min(PA(:));
max1 = max(PA(:));
min2 = min(NA(:));
max2 = max(NA(:));

disp(n_spikes_pos)
disp(n_spikes_neg)
disp(total_spikes)

for i = 1:30
    subplot(2,1,1)
    imagesc(PA(:,:,i)');
    axis image;
    colormap gray;
    caxis([min1 max1])
    title('Spike')
    subplot(2,1,2)
    imagesc(NA(:,:,i)');
    axis image;
    colormap gray;
    caxis([min2 max2])
    title('No Spike')
    pause(1)
end


end