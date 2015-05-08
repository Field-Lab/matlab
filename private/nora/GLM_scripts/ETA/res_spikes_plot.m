function res_spikes_plot(testmovie, res)

frames = length(testmovie);
writerObj = VideoWriter('/Volumes/Lab/Users/Nora/res_Spikes.avi');
open(writerObj);

for i_frame = 1:frames
    disp(i_frame)
    imagesc(testmovie(:,:,i_frame)')
    axis image
    colormap gray
    axis off
    hold on
    idx = (i_frame-1)*10+(1:10);
    response = sum(res.spikes(:,idx),2);
    scatter(res.centers(:,2), res.centers(:,1), 1+100*response, 'r');
    F = getframe;
    writeVideo(writerObj, F);
    clf
end

close(writerObj)

end
