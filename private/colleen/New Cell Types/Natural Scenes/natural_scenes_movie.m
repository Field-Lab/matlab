num_of_images = 20;
v = VideoWriter('NSEM movie');
v.FrameRate = 30; 
open(v)

for i = 1:num_of_images

movie = load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian/matfiles/movie_chunk_', num2str(i+849),'.mat']);

A = permute(movie.movie(:,:,1:4:end), [1 2 4 3])/255;
writeVideo(v,A);

end

close(v)