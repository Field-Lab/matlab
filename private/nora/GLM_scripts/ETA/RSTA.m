function [ETA] = RSTA(xval, movie, threshold)

n_blocks = length(movie);
stim_size = size(movie{1,1}.matrix);

n_spikes = 0;
ETA = zeros(stim_size(1), stim_size(2), 90);
%avg = zeros(stim_size(1), stim_size(2), 30);

for i_block = 1:n_blocks
    mov = movie{i_block}.matrix;
    res_spikes = find(logical(xval{i_block}.rasters.recorded) & (xval{i_block}.glm_rateperbin < threshold*mean(xval{i_block}.glm_rateperbin)));
 %   avg_stim = find(xval{i_block}.glm_rateperbin < threshold*mean(xval{i_block}.glm_rateperbin));
    
    n_spikes = n_spikes + length(res_spikes);
    for i_spike = 1:length(res_spikes)
        frame = ceil(res_spikes(i_spike)/10);
        if frame>89
            ETA = ETA + double(mov(:,:,(frame-89):frame));
        end
    end
%     for i_spike = 1:length(avg_stim)
%         frame = ceil(avg_stim(i_spike)/10);
%         if frame>29
%             avg = avg + double(mov(:,:,(frame-29):frame));
%         end
%     end
%     
end


%%
if 0
for i = 1:30
    subplot(1,2,1)
    imagesc(-ETA(:,:,i))
    axis image
    colormap gray
    caxis([5.6 7.7])
    title('Error Triggered Average, across Cells and Time')
    subplot(1,2,2)
    imagesc(avg_image(:,:,i))
    axis image
    colormap gray
    caxis([1.85 2.65])
    title('Average Image')
    pause(0.1)
end
end

disp(n_spikes)






    









ETA = ETA/n_spikes;
% 
% min1 = min(ETA(:));
% max1 = max(ETA(:));
% 
% % %%
% % h = figure(1);
% % 
% % % for i = 1:90
% % %     
% % %     imagesc(ETA(:,:,i)');
% % %     axis image;
% % %     colormap gray;
% % %     caxis([min1 max1])
% % %     title('Residual STA')
% % %     pause(0.1)
% % % 
% % % end
>>>>>>> 6daf92d2f0b7134b21a9d2e54a9f1d77f811b7a3:private/nora/GLM_scripts/ETA/RSTA.m

end

