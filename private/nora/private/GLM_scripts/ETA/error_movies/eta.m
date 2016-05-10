ROI = 5;
ETA = zeros(2*ROI+1, 2*ROI+1, 30);
avg_image = zeros(2*ROI+1, 2*ROI+1, 30);
edge = 0;
n_bins = length(res.spikes);

for i_cell = 1:length(res.cells)
    center = res.centers(i_cell);
    try
        relevant_stim = double(testmovie((res.centers(2)-ROI):(res.centers(2)+ROI), (res.centers(1)-ROI):(res.centers(1)+ROI), :));
        for i_bin = 1:n_bins
            frame = ceil(i_bin/10);
            if frame >29
                ETA = ETA + res.spikes(i_cell, i_bin)*relevant_stim(:,:,(frame-29):frame);
                avg_image = avg_image + relevant_stim(:,:,(frame-29):frame);
            end
        end
    catch
        edge = edge +1;
        warning(['cell ' num2str(i_cell) ' is too close to edge of stim'])
    end
    disp(i_cell)
end
save('/Users/Nora/Desktop/ETA11.mat', 'ETA')
save('/Users/Nora/Desktop/avg11.mat', 'avg_image')

%%
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