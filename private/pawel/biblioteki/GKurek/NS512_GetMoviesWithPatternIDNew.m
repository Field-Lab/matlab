function SelectedMovies = NS512_GetMoviesWithPatternIDNew(full_path,patternID);

SelectedMovies = [];

NS_GlobalConstants = NS_GenerateGlobalConstants(512);
header=readMHchunk_GlobalPath(full_path);
NumberOfMovies=header.number_of_movies;
for x=1:NumberOfMovies
    MovieData4=NS_MovieData_GlobalPath(full_path,x,NS_GlobalConstants);
    
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);
    patternIDList = MovieData(2:3:end,:);
    ind = find(patternIDList==patternID);
    
    if (size(ind,1)~=0) 
        SelectedMovies = [ SelectedMovies; x]; 
    end
end