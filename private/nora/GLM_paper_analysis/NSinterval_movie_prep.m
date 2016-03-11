movie_path = '/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/';
fitmovie = zeros(40, 80, 3600*120+500, 'uint8');
i_start = 1;
i_chunk = 1;
tic
while i_start < 3600*120;
    chunk = load([movie_path 'movie_chunk_' num2str(i_chunk) '.mat']);
    chunk_length = size(chunk.movie_chunk, 3);
if chunk_length>0    
    i_end = i_start + chunk_length -1;
    fitmovie(:,:,i_start:i_end) = permute(imresize(chunk.movie_chunk,1/4, 'box'), [2, 1, 3]);
    i_start = i_end+1;
end   
 i_chunk = i_chunk +1;
if ~mod(i_chunk, 10)
disp(i_end)
end
end
toc
clear chunk
fitmovie = fitmovie(:,:,1:3600*120);

save('/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/downsampledNSinterval.mat', fitmovie, '-v7.3');
