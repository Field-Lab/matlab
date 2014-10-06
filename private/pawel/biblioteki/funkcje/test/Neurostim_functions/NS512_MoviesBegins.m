function MoviesBegins=NS512_MoviesBegins(full_path,NS_GlobalConstants);

fid0=fopen(full_path,'r','b');
header=readMHchunk(fid0);
fclose(fid0);
NumberOfMovies=header.number_of_movies;
MoviesBegins=zeros(1,NumberOfMovies);

for i=1:NumberOfMovies
    MovieData=NS_MovieData_GlobalPath(full_path,i,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData);
    MoviesBegins(1,i)=MovieBegin;
end
