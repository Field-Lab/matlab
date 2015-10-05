function res_spikes_plot(testmovie, res, filename)

default_colors = get(gca,'ColorOrder');
frames = length(testmovie);
writerObj = VideoWriter(filename);
writerObj.FrameRate = 20;
open(writerObj);

n_cells = length(res.centers);
colors = zeros(n_cells,3);

for i_frame = 1:frames/20
    if ~mod(i_frame,100); disp(i_frame); end
    imagesc(testmovie(:,:,i_frame)')
    axis image
    colormap gray
    % axis off
    hold on
    % idx = (i_frame-1)*10+(1:10);
    idx = i_frame;
    response = sum(res.spikes(:,idx),2);
    
    neg = sum(response<0);
    pos = sum(response>0) ;
    
    colors(response<0,:) = repmat(default_colors(1,:), neg, 1);
    colors(response>0,:) = repmat(default_colors(2,:), pos, 1);
    
    scatter(res.centers(:,2), res.centers(:,1), 1+10*abs(response), colors, 'filled');
    
%     drawnow
%     ax=gca;
%     ax.Units = 'pixels';
%     ax.Position
%     
    F = getframe(gcf, [73.8000 108.2  434.0000  216.6]);
    writeVideo(writerObj, F);
    clf
end

close(writerObj)

end
