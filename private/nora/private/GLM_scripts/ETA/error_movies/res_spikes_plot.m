function res_spikes_plot(testmovie, res, filename, varargin)

p = inputParser;
p.addParameter('scaling', 1) % use this is the classification is not at the same res as the movie you want
p.addParameter('frame_rate', 20);
p.parse(varargin{:});

res.centers = res.centers * p.Results.scaling;

default_colors = get(gca,'ColorOrder');
frames = length(testmovie);
writerObj = VideoWriter(filename);
writerObj.FrameRate = p.Results.frame_rate;
open(writerObj);

n_cells = length(res.centers);
colors = zeros(n_cells,3);

largest = max(res.spikes(:));
dot_scaling = 500/largest;

for i_frame = 1:frames
    if ~mod(i_frame,100); disp(i_frame); end
    imagesc(testmovie(:,:,i_frame))
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
    
    scatter(res.centers(:,1), size(testmovie,1) - res.centers(:,2), 5+dot_scaling*abs(response), default_colors(1,:), 'filled');
    F = getframe(gcf, [73.8000 108.2  434.0000  216.6]);
    writeVideo(writerObj, F);
    clf
end

close(writerObj)

end
