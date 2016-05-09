function plot_movies(recon_movie, original_movie, varargin)

% both movies should be size (space1, space1, time)


p = inputParser;
p.addParameter('fraction_of_frames', 1, @isnumeric)
p.addParameter('pause_length', 'click')
p.parse(varargin{:});

movie_size = size(original_movie);

figure;

frames = min(movie_size(3),length(recon_movie))*p.Results.fraction_of_frames;
for i_frame = 1:frames
    
plot(original_movie(40,:,i_frame))
hold on
plot(recon_movie(40,:,i_frame))
hold off

if isnumeric(p.Results.pause_length)
    pause(p.Results.pause_length)
else
    pause()
end

end



end