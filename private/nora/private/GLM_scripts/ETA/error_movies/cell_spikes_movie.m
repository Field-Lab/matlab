function cell_spikes_movie(testmovie, res, filename, datarun)

default_colors = get(gca,'ColorOrder');
frames = length(testmovie);
% writerObj_stim = VideoWriter(filename);
writerObj_spike = VideoWriter(filename);
% writerObj_stim.FrameRate = 10;
writerObj_spike.FrameRate = 15;
% open(writerObj_stim);
open(writerObj_spike);

% n_cells = length(res.cid);
% colors = zeros(n_cells,3);

% for i_frame = 1:frames
%     if ~mod(i_frame,100); disp(i_frame); end
%     
%     figure(1);
%     imagesc(testmovie(:,:,i_frame)')
%     axis image
%     colormap gray
%     % axis off
%     % hold on
%     % idx = (i_frame-1)*10+(1:10);
%     % idx = i_frame;
%     % response = sum(res.spikes(:,idx),2);
%     
%     % neg = sum(response<0);
%     %  pos = sum(response>0) ;
%     
%     % colors(response<0,:) = repmat(default_colors(1,:), neg, 1);
%     % colors(response>0,:) = repmat(default_colors(2,:), pos, 1);
%     F = getframe(gcf, [73.8000 108.2  434.0000  216.6]);
%     writeVideo(writerObj_stim, F);
%     clf
%     
% end
for i_frame = 1:frames
    if ~mod(i_frame,100); disp(i_frame); end
    
    plot_rf_fit(datarun, res.cid, 'fill', true, 'alpha', 0.5, 'fill_color', default_colors(1,:), 'edge', false)
    spiking_cells = res.cid(logical(res.spikes(:,i_frame)));
    hold on
    plot_rf_fit(datarun, spiking_cells, 'fill', true, 'alpha', 0.5, 'fill_color', default_colors(2,:), 'edge', false)
    hold off

%    scatter(res.centers(:,2), res.centers(:,1), 1+10*abs(response), colors, 'filled');
    
%     drawnow


    axis_locations = [73.8000 80.2 434.0000  246.6];
    ax=gca;
    ax.Units = 'pixels';
    ax.Position = axis_locations;
    axis off
%     
    F = getframe(gcf, axis_locations);
    writeVideo(writerObj_spike, F);
    clf
end

close(writerObj_spike)
%close(writerObj_stim)

end
