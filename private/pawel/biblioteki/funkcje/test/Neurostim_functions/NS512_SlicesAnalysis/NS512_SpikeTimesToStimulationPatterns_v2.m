function [MovieNumber,RepetitionNumber,PatternNumber,Latency,PatternID]=NS512_SpikeTimesToStimulationPatterns_v2(full_path,SpikeTime,MoviesBegins,NS_GlobalConstants);

a=find((SpikeTime-MoviesBegins)>0);
if length(a)==0
    MovieNumber=0;
    RepetitionNumber=0;
    PatternNumber=0;
    Latency=0;
    PatternID=0;
    return;
else
    MovieID=a(length(a));
    MovieNumber=MovieID;
end

TimeRelativeToMovieBegin=SpikeTime-MoviesBegins(MovieID);
MovieDataFull=NS_MovieData_GlobalPath(full_path,MovieID,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieDataFull);

RepetitionNumber=ceil(TimeRelativeToMovieBegin/RepetPeriod);
if RepetitionNumber>RepetNumber
    MovieNumber=0;
    RepetitionNumber=0;
    PatternNumber=0;
    Latency=0;
    PatternID=0;
    return;
end

TimeRelativeToRepetitionBegin=TimeRelativeToMovieBegin-(RepetitionNumber-1)*RepetPeriod;
patternTimeList = MovieData(1:3:end,:);
PatternsList=MovieData(2:3:end,:);

a=find((TimeRelativeToRepetitionBegin-patternTimeList)>0);
if length(a)==0
    MovieNumber=0;
    RepetitionNumber=0;
    PatternNumber=0;
    Latency=0;
    PatternID=0;
    return;
else
    PatternID=a(length(a));
    PatternNumber=MovieData(2+(PatternID-1)*3);
    %AllAppearanceOfThiPattern=find(
    Latency=TimeRelativeToRepetitionBegin-MovieData(1+(PatternID-1)*3);
end