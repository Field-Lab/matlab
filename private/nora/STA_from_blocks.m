function STA = STA_from_blocks(fitspikes, fitmovie)
    stim_size = size(fitmovie{1});
    STA = zeros(stim_size(1), stim_size(2), 30); % 30 frame STA
    for i_block = 1:length(fitmovie)
        block_spikes = fitspikes{i_block, 1};
        for i_spike = 1:length(block_spikes)
            spike_frame = floor(120*block_spikes(i_spike));
            if spike_frame > 29
                STA = STA + double(fitmovie{i_block}(:,:,(spike_frame-29):spike_frame));
            end
        end
    end
    figure;
    for i = 1:30
        imagesc(STA(:,:,i))
        axis image
        colormap gray
        pause(0.1)
    end
    imagesc(sum(STA(:,:,24:27), 3))
    axis image
    colormap gray
end