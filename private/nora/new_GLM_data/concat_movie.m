function fitmovie = concat_movie(fitmovie_cell)
n_blocks = length(fitmovie_cell);
fitmovie = [];
for block = 1:n_blocks
    fitmovie = cat(3, fitmovie, fitmovie_cell{block});
end
end