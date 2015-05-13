
min1 = min(eta(:));
max1 = max(eta(:));
min2 = min(sta(:));
max2 = max(sta(:));

writerObj = VideoWriter('/Users/Nora/Desktop/resSTA_ON2.avi');
writerObj.FrameRate = 20;
open(writerObj);

h = figure(1);

for i = 1:90
    subplot(2,1,1)
    imagesc(eta(:,:,i)');
    axis image;
    colormap gray;
    caxis([min1 max1])
    title('Residual STA')
    subplot(2,1,2)
    imagesc(sta(:,:,i)');
    axis image;
    colormap gray;
    caxis([min2 max2])
    title('STA')
    F = getframe(h);
    writeVideo(writerObj, F);
    clf
end

close(writerObj)