for i_STA = 1:30
    subplot(1,2,1)
    imagesc(squeeze(STAorig2(i_STA,ROI.xvals,ROI.yvals)))
    colormap gray
    axis image
    subplot(1,2,2)
    imagesc(squeeze(STA_concat(:,:,i_STA)))
    colormap gray
    axis image
    pause(0.1)
end