ROI = 11;
ETA = zeros(2*ROI+1, 2*ROI+1, 90);
avg_image = zeros(2*ROI+1, 2*ROI+1, 90);
edge = 0;
n_bins = length(res.spikes);

for i_cell = 1:length(res.cells)
    center = res.centers(i_cell);
    try
        relevant_stim = double(testmovie((res.centers(2)-ROI):(res.centers(2)+ROI), (res.centers(1)-ROI):(res.centers(1)+ROI), :));
        for i_bin = 1:n_bins
            frame = ceil(i_bin/10);
if frame > 89
                ETA = ETA + res.spikes(i_cell, i_bin)*relevant_stim(:,:,(frame-89):frame);
                avg_image = avg_image + relevant_stim(:,:,(frame-89):frame);
            end
        end
    catch
        edge = edge +1;
        warning(['cell ' num2str(i_cell) ' is too close to edge of stim'])
    end
    disp(i_cell)
end
save('/Volumes/Lab/Users/Nora/OFFETA11_90.mat', 'ETA')
save('/Users/Nora/Desktop/OFFavg11_90.mat', 'avg_image')

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
