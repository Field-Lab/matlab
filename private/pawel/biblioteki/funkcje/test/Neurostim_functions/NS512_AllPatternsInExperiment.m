function AllPatternsFinal=NS512_AllPatternsInExperiment(full_path,NS_GlobalConstants);

fid0=fopen(full_path,'r','b');
header=readMHchunk(fid0);
fclose(fid0);
NumberOfMovies=header.number_of_movies;
MoviesBegins=zeros(1,NumberOfMovies);

AllPatterns=[];
for i=1:NumberOfMovies
    MovieDataAll=NS_MovieData_GlobalPath(full_path,i,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieDataAll);
    PatternsInMovie=MovieData(2:3:length(MovieData))
    AllPatterns=[AllPatterns' PatternsInMovie']';
    %MoviesBegins(1,i)=MovieBegin;
end
AllPatternsFinal=unique(AllPatterns);