movie_path = '/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian/matfiles/';
tic
mov_letter = {'A', 'B', 'C', 'D', 'E', 'F'};
offset = 0;
for movie = 1:6
    fitmovie = zeros(40, 80, 3000*120, 'uint8');
    idx = 1:120;
    for i_chunk = (1:3000)+offset
        chunk = load([movie_path 'movie_chunk_' num2str(i_chunk) '.mat']);
        fitmovie(:,:,idx) = permute(imresize(chunk.movie,1/4, 'box'), [2, 1, 3]);
        idx = idx+120;
        if ~mod(i_chunk, 100)
            disp(i_chunk)
        end
    end
    toc
    save(['/Volumes/Lab/Users/Nora/downsampledNSbrownian_' mov_letter{movie} '_3000.mat'], 'fitmovie', '-v7.3');
end