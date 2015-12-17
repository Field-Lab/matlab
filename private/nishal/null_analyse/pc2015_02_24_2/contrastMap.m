function cMap=contrastMap(condMovies)
condMovies = permute(condMovies,[2,3,1]);
cMap  = sqrt(mean(condMovies.^2,3));
end