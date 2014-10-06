function [NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,Patterns]=NS512_MoviePatterns(MovieFilePath,NS_GlobalConstants);
% NumberOfMovies - how many movies exists in given movie file
% NumberOfPatternsPerMovie - how many different patters are used in each
% movie. This is vector with length equal to NumberOfMovies.
% AllPatternsUsed - all patterns used in ALL movies. For number of
% different patterns used in the experiment, just type
% length(AllPatternsUsed).
% Patterns - NxM array, where N=NumberOfMovies, M=512. This gives
% information how many times given pattern was used in given movie.

fid0=fopen(MovieFilePath,'r','b');
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies;
fclose(fid0);

Patterns=zeros(NumberOfMovies,512);
for i=1:NumberOfMovies
    MovieData1=NS512_MovieData2(MovieFilePath,i,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData1);
    Pt=MovieData(2:3:length(MovieData)); %all patterns in given movie
    for j=1:length(Pt)
        Patterns(i,Pt(j))=Patterns(i,Pt(j))+1;
    end
end
NumberOfPatternsPerMovie=sum(Patterns,2);
AllPatternsUsed=find(sum(Patterns));