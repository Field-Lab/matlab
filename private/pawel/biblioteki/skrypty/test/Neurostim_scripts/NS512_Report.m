ChipAddresses=[24:31];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=500;

cd G:\2010-09-11-0;

TotalNumberOfPatterns=512;
TotalNumberOfMovies=224;

PatternsForMovies=zeros(TotalNumberOfMovies,TotalNumberOfPatterns);
for i=1:TotalNumberOfMovies
    ChunkData=NS_MovieData('005',i,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
    Patterns=MovieData(2:3:length(MovieData));
    PatternsForMovies(i,Patterns)=1; 
    %plot(ones(1,length(Patterns))*i,Patterns,'bd')
    %hold on
end
    
    