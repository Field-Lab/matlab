function display_movies(recon_movie, original_movie, varargin)

% both movies should be size (space1, space1, time)


p = inputParser;
p.addParameter('fraction_of_frames', 1, @isnumeric)
p.addParameter('pause_length', 0.1, @isnumeric)
p.addParameter('disp_range', 0)
p.parse(varargin{:});

movie_size = size(original_movie);

if any(p.Results.disp_range)
    disp_range = p.Results.disp_range;
else
    disp_range = [0 movie_size(1) 0 movie_size(2)];
end

figure;

frames = min(movie_size(3),length(recon_movie))*p.Results.fraction_of_frames;
for i_frame = 1:frames
    
subplot(2,1,1)
imagesc(original_movie(:,:,i_frame)')
colormap gray
caxis([-0.5 0.5])
axis image
axis(disp_range)
axis off


subplot(2,1,2)
imagesc(recon_movie(:,:,i_frame)')
colormap gray
caxis([-0.5 0.5])
axis image
axis(disp_range)
axis off

pause(p.Results.pause_length)

end



end