function [X]=NS512_MovieNumberAndTimingForPattern(FilePath,PatternNumber);
NS_GlobalConstants=NS_GenerateGlobalConstants(512);

[NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,Patterns]=NS512_MoviePatterns(FilePath,NS_GlobalConstants);
Movies=find(Patterns(:,PatternNumber)>0);

X=zeros(length(Movies),3);
X(:,1)=Movies;
for i=1:length(Movies)
    ChunkData=NS_MovieData_GlobalPath(FilePath,Movies(i),NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
    index=find(MovieData(2:3:length(MovieData))==PatternNumber);
    X(i,2)=length(index);
    for j=1:length(index)
        Time=MovieData(index(j)*3-2);
        X(i,j+2)=Time;
    end
end